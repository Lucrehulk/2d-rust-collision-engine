// Room size.
const ROOM_SIZE: f32 = 1024.0;
// Amount of spatial grid units on a single side (i.e. area = SPATIAL_GRID_DIMENSION * SPATIAL_GRID_DIMENSION).
const SPATIAL_GRID_DIMENSION: usize = 64;
// Target number of threads for the application to use.
const THREADS: usize = 16;
// The capacity for the maximum number of entities the system is designed to handle. If this value is exceeded, the engine will fail.
const MAX_ENTITIES: usize = 10000;
// The capacity for the maximum number of open entity indexes the system is designed to handle. if this value is exceeded, the engine will fail.
const MAX_ENTITIES_TO_REPLACE: usize = 50;
// Determines if friction is enabled.
const FRICTION: bool = false;
// The constant of friction, this is consistently multiplied to the entity's x and y velocity every tick.
const FRICTIONAL_CONSTANT: f32 = 0.99;
// Determines if gravity is enabled.
const GRAVITY: bool = false;
// The constant of the acceleration due to gravity, this is consistently added to the entity's y velocity so as to make the entity go downwards (hence replicating gravity).
const GRAVITY_CONSTANT: f32 = 0.2;
// The acceleration applied to both the x and y velocities when a collision occurs. 
const COLLISION_ACCELERATION_CONSTANT: f32 = 0.1;
// The acceleration applied to an entity when it is in a movement state.
const MOVEMENT_ACCELERATION_CONSTANT: f32 = 0.01;
// The downtime, in ms, per tick.
const TICK_TIME: u64 = 16;

// DO NOT SET THESE
const ENCODING_BITS: usize = (usize::BITS - SPATIAL_GRID_DIMENSION.leading_zeros() - 1) as usize;
const SPATIAL_GRID_AREA: usize = SPATIAL_GRID_DIMENSION * SPATIAL_GRID_DIMENSION;
const ROOM_GRID_RATIO: f32 = (SPATIAL_GRID_DIMENSION as f32) / ROOM_SIZE;
const MOVEMENT_ACCELERATION_CONSTANT_45_DEG: f32 = MOVEMENT_ACCELERATION_CONSTANT * 1.4142135 / 2.0;

use std::ptr;
use std::collections::HashSet;
use std::sync::Mutex;

struct Entity {
    index: usize,
    grid_pos_x: usize,
    grid_pos_y: usize,
    grid_body: usize,
    x: f32,
    y: f32,
    velocity_x: f32,
    velocity_y: f32,
    max_velocity_x: f32,
    max_velocity_y: f32,
    acceleration_x: f32,
    acceleration_y: f32,
    radius: f32,
    body_type: u8,
    movable: bool,
    replace: bool
}

struct Room {
    spatial_grid: Vec<Mutex<HashSet<usize>>>,
    spatial_grid_ptr: *mut Mutex<HashSet<usize>>,
    collision_positions: Mutex<HashSet<usize>>,
    collision_positions_ptr: *mut Mutex<HashSet<usize>>,
    entities: Vec<Entity>,
    entities_ptr: *mut Entity,
    chunks: [[usize; 2]; THREADS],
    chunks_count: usize,
    replacement_queue: Vec<usize>,
    replacement_queue_ptr: *mut usize,
    tick: usize
}

#[inline(always)]
fn manage_collision(dx: f32, dy: f32, dsq: f32, entity: &mut Entity, collision_entity: &mut Entity) {
    let inverse_distance = 1.0 / dsq.sqrt();
    let nx = dx * inverse_distance;
    let ny = dy * inverse_distance;
    entity.velocity_x -= nx * COLLISION_ACCELERATION_CONSTANT;
    entity.velocity_y -= ny * COLLISION_ACCELERATION_CONSTANT;
    collision_entity.velocity_x += nx * COLLISION_ACCELERATION_CONSTANT;
    collision_entity.velocity_y += ny * COLLISION_ACCELERATION_CONSTANT;
}

#[inline(always)]
fn update_grid_position(pos: usize, index: usize, spatial_grid: *mut Mutex<HashSet<usize>>, collision_positions: *mut Mutex<HashSet<usize>>) {
    unsafe { 
        let mut position = (*spatial_grid.add(pos)).lock().unwrap();
        if position.len() == 1 {
            (*collision_positions).lock().unwrap().insert(pos); 
        };
        position.insert(index); 
    };
}

#[inline(always)]
fn remove_grid_position(pos: usize, index: &usize, spatial_grid: *mut Mutex<HashSet<usize>>, collision_positions: *mut Mutex<HashSet<usize>>) {
    unsafe { 
        let mut position = (*spatial_grid.add(pos)).lock().unwrap();
        if position.len() == 2 {
            (*collision_positions).lock().unwrap().remove(&pos);
        };
        position.remove(index); 
    };
}

#[inline(always)]
fn update_entity_body(y_bound: usize, x_bound: usize, y_position_update: usize, x_position_update: usize, y_position_delete: usize, x_position_delete: usize, index: usize, spatial_grid: *mut Mutex<HashSet<usize>>, collision_positions: *mut Mutex<HashSet<usize>>) {
    for y_offset in y_position_update..usize::min(y_position_update + y_bound, SPATIAL_GRID_DIMENSION) {
        for x_offset in x_position_update..usize::min(x_position_update + x_bound, SPATIAL_GRID_DIMENSION) {
            let position = (y_offset << ENCODING_BITS) | x_offset;
            update_grid_position(position, index, spatial_grid, collision_positions);
        };
    };
    for y_offset in y_position_delete..usize::min(y_position_delete + y_bound, SPATIAL_GRID_DIMENSION) {
        for x_offset in x_position_delete..usize::min(x_position_delete + x_bound, SPATIAL_GRID_DIMENSION) {
            let dposition = (y_offset << ENCODING_BITS) | x_offset;
            remove_grid_position(dposition, &index, spatial_grid, collision_positions);
        };
    };
}

impl Room {
    fn init() -> Room {
        let mut entities = Vec::with_capacity(MAX_ENTITIES);
        let entities_ptr = entities.as_mut_ptr();
        let mut replacement_queue = Vec::with_capacity(MAX_ENTITIES_TO_REPLACE);
        let replacement_queue_ptr = replacement_queue.as_mut_ptr();
        let mut spatial_grid = Vec::with_capacity(SPATIAL_GRID_AREA);
        for _pos in 0..SPATIAL_GRID_AREA {
            spatial_grid.push(Mutex::new(HashSet::new()));
        };
        let spatial_grid_ptr = spatial_grid.as_mut_ptr();
        let collision_positions = Mutex::new(HashSet::with_capacity(SPATIAL_GRID_AREA));
        Room {
            spatial_grid,
            spatial_grid_ptr,
            collision_positions,
            collision_positions_ptr: ptr::null_mut(),
            entities,
            entities_ptr,
            chunks: [[0, 0]; THREADS],
            chunks_count: 0,
            replacement_queue,
            replacement_queue_ptr,
            tick: 0
        }
    }

    fn update_chunks(&mut self) {
        let len = self.entities.len();
        if len > THREADS {
            let chunk_size = len / THREADS;
            let chunk_rem = len % THREADS;
            let mut chunk_pos = 0;
            for chunk in 0..THREADS {
                let next_chunk = if chunk < chunk_rem {
                    chunk_pos + chunk_size + 1
                } else {
                    chunk_pos + chunk_size
                };
                self.chunks[chunk] = [chunk_pos, next_chunk];
                chunk_pos = next_chunk;
            };
            self.chunks_count = THREADS;
        } else {
            for chunk in 0..len {
                self.chunks[chunk] = [chunk, chunk + 1];
            };
            self.chunks_count = len;
        };
    }

    fn create_entity(&mut self, x: f32, y: f32, velocity_x: f32, velocity_y: f32, max_velocity_x: f32, max_velocity_y: f32, radius: f32, body_type: u8) {
        let grid_body = (radius * 2.0 * ROOM_GRID_RATIO).ceil() as usize + 1;
        let grid_pos_x = ((x - radius) * ROOM_GRID_RATIO) as usize;
        let grid_pos_y = ((y - radius) * ROOM_GRID_RATIO) as usize;
        let mut index = usize::MAX;
        if self.replacement_queue.len() == 0 {
            unsafe {
                index = self.entities.len();
                *self.entities_ptr.add(self.entities.len()) = Entity {
                    index,
                    grid_pos_x,
                    grid_pos_y,
                    grid_body,
                    x,
                    y,
                    velocity_x,
                    velocity_y,
                    max_velocity_x,
                    max_velocity_y,
                    acceleration_x: 0.0,
                    acceleration_y: 0.0,
                    radius,
                    body_type,
                    movable: max_velocity_x != 0.0 || max_velocity_y != 0.0,
                    replace: false
                };
                self.entities.set_len(self.entities.len() + 1);
            };
        } else {
            unsafe { 
                index = *self.replacement_queue_ptr;
                *self.entities_ptr.add(index) = Entity {
                    index,
                    grid_pos_x,
                    grid_pos_y,
                    grid_body,
                    x,
                    y,
                    velocity_x,
                    velocity_y,
                    max_velocity_x,
                    max_velocity_y,
                    acceleration_x: 0.0,
                    acceleration_y: 0.0,
                    radius,
                    body_type,
                    movable: max_velocity_x != 0.0 || max_velocity_y != 0.0,
                    replace: false
                }; 
                let len = self.replacement_queue.len() - 1;
                ptr::swap(self.replacement_queue_ptr, self.replacement_queue_ptr.add(len));
                self.replacement_queue.set_len(len);
            };
        };
        for y_pos in grid_pos_y..usize::min(grid_pos_y + grid_body, SPATIAL_GRID_DIMENSION) {
            for x_pos in grid_pos_x..usize::min(grid_pos_x + grid_body, SPATIAL_GRID_DIMENSION) {
                update_grid_position((y_pos << ENCODING_BITS) | x_pos, index, self.spatial_grid_ptr, self.collision_positions_ptr);
            };
        };
        self.update_chunks();
    }

    fn remove_entity(&mut self, index: usize) {
        unsafe { 
            let entity = self.entities.get_unchecked_mut(index);
            entity.replace = true;
            *self.replacement_queue_ptr.add(self.replacement_queue.len()) = entity.index;
            self.replacement_queue.set_len(self.replacement_queue.len() + 1);
            for y_pos in entity.grid_pos_y..usize::min(entity.grid_pos_y + entity.grid_body, SPATIAL_GRID_DIMENSION) {
                for x_pos in entity.grid_pos_x..usize::min(entity.grid_pos_x + entity.grid_body, SPATIAL_GRID_DIMENSION) {
                    remove_grid_position((y_pos << ENCODING_BITS) | x_pos, &entity.index, self.spatial_grid_ptr, self.collision_positions_ptr);
                };
            };
        };
        self.update_chunks();
    }

    fn create_entity_movement_from_angle(&mut self, index: usize, angle: f32) {
        let entity = unsafe { self.entities.get_unchecked_mut(index) };
        entity.acceleration_x = f32::cos(angle) * MOVEMENT_ACCELERATION_CONSTANT;
        entity.acceleration_y = f32::sin(angle) * MOVEMENT_ACCELERATION_CONSTANT;
    }

    fn create_entity_movement_from_cardinal_direction(&mut self, index: usize, direction: u8) {
        let entity = unsafe { self.entities.get_unchecked_mut(index) };
        match direction {
            0 => {
                entity.acceleration_x = MOVEMENT_ACCELERATION_CONSTANT;
            }
            1 => {
                entity.acceleration_x = -MOVEMENT_ACCELERATION_CONSTANT;
            }
            2 => {
                entity.acceleration_y = MOVEMENT_ACCELERATION_CONSTANT;
            }
            3 => {
                entity.acceleration_y = -MOVEMENT_ACCELERATION_CONSTANT;
            }
            4 => {
                entity.acceleration_x = MOVEMENT_ACCELERATION_CONSTANT_45_DEG;
                entity.acceleration_y = MOVEMENT_ACCELERATION_CONSTANT_45_DEG;
            }
            5 => {
                entity.acceleration_x = MOVEMENT_ACCELERATION_CONSTANT_45_DEG;
                entity.acceleration_y = -MOVEMENT_ACCELERATION_CONSTANT_45_DEG;
            }
            6 => {
                entity.acceleration_x = -MOVEMENT_ACCELERATION_CONSTANT_45_DEG;
                entity.acceleration_y = MOVEMENT_ACCELERATION_CONSTANT_45_DEG;
            }
            7 => {
                entity.acceleration_x = -MOVEMENT_ACCELERATION_CONSTANT_45_DEG;
                entity.acceleration_y = -MOVEMENT_ACCELERATION_CONSTANT_45_DEG;
            }
            _ => {}
        }
    }

    fn stop_entity_movement(&mut self, index: usize) {
        let entity = unsafe { self.entities.get_unchecked_mut(index) };
        entity.acceleration_x = 0.0;
        entity.acceleration_y = 0.0;
    }

    fn stop_entity_movement_x(&mut self, index: usize) {
        let entity = unsafe { self.entities.get_unchecked_mut(index) };
        entity.acceleration_x = 0.0;
    }

    fn stop_entity_movement_y(&mut self, index: usize) {
        let entity = unsafe { self.entities.get_unchecked_mut(index) };
        entity.acceleration_y = 0.0;
    }

    fn update(&mut self, entities_ptr: usize, spatial_grid_ptr: usize, collision_positions_ptr: usize) {
        let mut spatial_grid_update_threads = Vec::with_capacity(self.chunks_count);
        for thread in 0..self.chunks_count {
            let chunk = self.chunks[thread];
            let worker = thread::spawn(move || {
                let entities = entities_ptr as *mut Entity;
                let spatial_grid = spatial_grid_ptr as *mut Mutex<HashSet<usize>>;
                let collision_positions = collision_positions_ptr as *mut Mutex<HashSet<usize>>;
                for index in chunk[0]..chunk[1] {
                    let entity = unsafe { &mut *entities.add(index) };
                    if !entity.replace {
                        if entity.movable {
                            if entity.acceleration_x != 0.0 {
                                entity.velocity_x += entity.acceleration_x;
                            };
                            if entity.acceleration_y != 0.0 {
                                entity.velocity_y += entity.acceleration_y;
                            };
                            if GRAVITY {
                                entity.velocity_y += GRAVITY_CONSTANT;
                            };
                            if FRICTION {
                                entity.velocity_x *= FRICTIONAL_CONSTANT;
                                entity.velocity_y *= FRICTIONAL_CONSTANT; 
                            };
                            if entity.velocity_x > entity.max_velocity_x {
                                entity.velocity_x = entity.max_velocity_x; 
                            } else if entity.velocity_x < -entity.max_velocity_x {
                                entity.velocity_x = -entity.max_velocity_x;
                            };
                            if entity.velocity_y > entity.max_velocity_y {
                                entity.velocity_y = entity.max_velocity_y; 
                            } else if entity.velocity_y < -entity.max_velocity_y {
                                entity.velocity_y = -entity.max_velocity_y;
                            };
                            entity.x += entity.velocity_x;
                            entity.y += entity.velocity_y; 
                            if entity.x - entity.radius < 0.0 {
                                entity.x = entity.radius;
                                entity.velocity_x = -entity.velocity_x;
                            };
                            if entity.y - entity.radius < 0.0 {
                                entity.y = entity.radius;
                                entity.velocity_y = -entity.velocity_y;
                            };
                            if entity.x + entity.radius > ROOM_SIZE {
                                entity.x = ROOM_SIZE - entity.radius;
                                entity.velocity_x = -entity.velocity_x;
                            };
                            if entity.y + entity.radius > ROOM_SIZE {
                                entity.y = ROOM_SIZE - entity.radius;
                                entity.velocity_y = -entity.velocity_y;
                            }; 
                            let spatial_grid_pos_x = ((entity.x - entity.radius) * ROOM_GRID_RATIO) as usize;
                            let spatial_grid_pos_y = ((entity.y - entity.radius) * ROOM_GRID_RATIO) as usize;
                            if spatial_grid_pos_x != entity.grid_pos_x || spatial_grid_pos_y != entity.grid_pos_y {
                                let shift_x = if spatial_grid_pos_x > entity.grid_pos_x {
                                    (spatial_grid_pos_x - entity.grid_pos_x, true)
                                } else {
                                    (entity.grid_pos_x - spatial_grid_pos_x, false)
                                };
                                let shift_y = if spatial_grid_pos_y > entity.grid_pos_y {
                                    (spatial_grid_pos_y - entity.grid_pos_y, true)
                                } else {
                                    (entity.grid_pos_y - spatial_grid_pos_y, false)
                                };
                                if shift_x.0 > entity.grid_body || shift_y.0 > entity.grid_body {
                                    update_entity_body(entity.grid_body, entity.grid_body, spatial_grid_pos_y, spatial_grid_pos_x, entity.grid_pos_y, entity.grid_pos_x, entity.index, spatial_grid, collision_positions);
                                } else {
                                    if shift_x.0 == 0 {
                                        if shift_y.1 {
                                            update_entity_body(shift_y.0, entity.grid_body, entity.grid_pos_y + entity.grid_body, spatial_grid_pos_x, entity.grid_pos_y, spatial_grid_pos_x, entity.index, spatial_grid, collision_positions);
                                        } else {
                                            update_entity_body(shift_y.0, entity.grid_body, spatial_grid_pos_y, spatial_grid_pos_x, spatial_grid_pos_y + entity.grid_body, spatial_grid_pos_x, entity.index, spatial_grid, collision_positions);
                                        };
                                    } else if shift_y.0 == 0 {
                                        if shift_x.1 {
                                            update_entity_body(entity.grid_body, shift_x.0, spatial_grid_pos_y, entity.grid_pos_x + entity.grid_body, spatial_grid_pos_y, entity.grid_pos_x, entity.index, spatial_grid, collision_positions);
                                        } else {
                                            update_entity_body(entity.grid_body, shift_x.0, spatial_grid_pos_y, spatial_grid_pos_x, spatial_grid_pos_y, spatial_grid_pos_x + entity.grid_body, entity.index, spatial_grid, collision_positions);
                                        };
                                    } else {
                                        if shift_y.1 {
                                            update_entity_body(shift_y.0, entity.grid_body, entity.grid_pos_y + entity.grid_body, spatial_grid_pos_x, entity.grid_pos_y, entity.grid_pos_x, entity.index, spatial_grid, collision_positions);
                                        } else {
                                            update_entity_body(shift_y.0, entity.grid_body, spatial_grid_pos_y, spatial_grid_pos_x, spatial_grid_pos_y + entity.grid_body, entity.grid_pos_x, entity.index, spatial_grid, collision_positions);
                                        };
                                        if shift_x.1 {
                                            update_entity_body(entity.grid_body - shift_y.0, shift_x.0, spatial_grid_pos_y, entity.grid_pos_x + entity.grid_body, spatial_grid_pos_y, entity.grid_pos_x, entity.index, spatial_grid, collision_positions);
                                        } else {
                                            update_entity_body(entity.grid_body - shift_y.0, shift_x.0, spatial_grid_pos_y, spatial_grid_pos_x, spatial_grid_pos_y,spatial_grid_pos_x + entity.grid_body, entity.index, spatial_grid, collision_positions);
                                        };
                                    };
                                };
                                entity.grid_pos_x = spatial_grid_pos_x;
                                entity.grid_pos_y = spatial_grid_pos_y; 
                            };
                        };
                    };
                };
            });
            spatial_grid_update_threads.push(worker);
        };
        for thread in spatial_grid_update_threads {
            thread.join().unwrap();
        };
        let mut collision_chunks = [[0, 0]; THREADS];
        let collision_positions = self.collision_positions.lock().unwrap();
        let collision_positions_vec = collision_positions.iter().collect::<Vec<&usize>>();
        let collision_positions_vec_ptr = collision_positions_vec.as_ptr() as usize;
        let mut len = collision_positions.len();
        if len > THREADS {
            let chunk_size = len / THREADS;
            let chunk_rem = len % THREADS;
            let mut chunk_pos = 0;
            for chunk in 0..THREADS {
                let next_chunk = if chunk < chunk_rem {
                    chunk_pos + chunk_size + 1
                } else {
                    chunk_pos + chunk_size
                };
                collision_chunks[chunk] = [chunk_pos, next_chunk];
                chunk_pos = next_chunk;
            };
            len = THREADS;
        } else {
            for chunk in 0..len {
                collision_chunks[chunk] = [chunk, chunk + 1];
            };
        };
        let mut collision_threads = Vec::with_capacity(len);
        for thread in 0..len {
            let chunk = collision_chunks[thread];
            let worker = thread::spawn(move || {
                let entities = entities_ptr as *mut Entity;
                let spatial_grid = spatial_grid_ptr as *mut Mutex<HashSet<usize>>;
                let collision_positions_vec = collision_positions_vec_ptr as *const &usize;
                for index in chunk[0]..chunk[1] {
                    let position = unsafe { **collision_positions_vec.add(index) }; 
                    let spatial_grid_position = unsafe { (*spatial_grid.add(position)).lock().unwrap() };
                    let mut entities_vec: Vec<*mut Entity> = Vec::with_capacity(spatial_grid_position.len());
                    for index in spatial_grid_position.iter().copied() {
                        let entity = unsafe { &mut *entities.add(index) };
                        for collision_index in &mut entities_vec {
                            let collision_entity = unsafe { &mut **collision_index };
                            if entity.body_type + collision_entity.body_type == 2 {
                                let dx = collision_entity.x - entity.x;
                                let dy = collision_entity.y - entity.y;
                                let dsq = dx * dx + dy * dy;
                                if dsq < (entity.radius + collision_entity.radius) * (entity.radius + collision_entity.radius) {
                                    manage_collision(dx, dy, dsq, entity, collision_entity);
                                };
                            } else if entity.body_type + collision_entity.body_type == 1 {
                                let (closest_x, closest_y, circle_center_x, circle_center_y, circle_radius) = if entity.body_type == 0 {
                                    let cx = collision_entity.x.clamp(entity.x - entity.radius, entity.x + entity.radius);
                                    let cy = collision_entity.y.clamp(entity.y - entity.radius, entity.y + entity.radius);
                                    (cx, cy, collision_entity.x, collision_entity.y, collision_entity.radius)
                                } else {
                                    let cx = entity.x.clamp(collision_entity.x - collision_entity.radius, collision_entity.x + collision_entity.radius);
                                    let cy = entity.y.clamp(collision_entity.y - collision_entity.radius, collision_entity.y + collision_entity.radius);
                                    (cx, cy, entity.x, entity.y, entity.radius)
                                };
                                if (closest_x - circle_center_x) * (closest_x - circle_center_x) + (closest_y - circle_center_y) * (closest_y - circle_center_y) < circle_radius * circle_radius {
                                    let dx = collision_entity.x - entity.x;
                                    let dy = collision_entity.y - entity.y;
                                    manage_collision(dx, dy, dx * dx + dy * dy, entity, collision_entity);
                                };
                            } else {
                                let dx = collision_entity.x - entity.x;
                                let dy = collision_entity.y - entity.y;
                                if f32::abs(dx) < entity.radius + collision_entity.radius && f32::abs(dy) < entity.radius + collision_entity.radius {
                                    manage_collision(dx, dy, dx * dx + dy * dy, entity, collision_entity);
                                };
                            };
                        };
                        entities_vec.push(entity as *mut Entity);
                    };
                };
            });
            collision_threads.push(worker);
        };
        for thread in collision_threads {
            thread.join().unwrap();
        };
        self.tick += 1;
    }
}

use std::{thread, time::{Duration, Instant}};

// Just example usage and stuff
fn main() {
    let mut world = Room::init();
    world.collision_positions_ptr = &mut world.collision_positions as *mut Mutex<HashSet<usize>>;
    let entities_ptr = world.entities_ptr as usize;
    let spatial_grid_ptr = world.spatial_grid_ptr as usize;
    let collision_positions_ptr = world.collision_positions_ptr as usize;
    
    // Create random circle entities example
    use rand::Rng; 
    let mut rng = rand::rng();

    for _ in 0..MAX_ENTITIES {
        let x = rng.random_range(0.0..ROOM_SIZE);
        let y = rng.random_range(0.0..ROOM_SIZE);
        let vx = rng.random_range(-2.0..2.0); 
        let vy = rng.random_range(-2.0..2.0); 
        let radius = rng.random_range(2.0..6.0);
        let body_type = /*f32::round(rng.random_range(0.0..1.0)) as u8*/ 1;
        world.create_entity(x, y, vx, vy, 2.0, 2.0, radius, body_type);
    }
    
    loop {
        let tick = Instant::now();
        world.update(entities_ptr, spatial_grid_ptr, collision_positions_ptr);
        let tick_time = tick.elapsed();
        println!("{:?}", tick_time);
        thread::sleep(Duration::from_millis(TICK_TIME));
    }
}
