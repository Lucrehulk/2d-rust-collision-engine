use super::config::config_data::{
    ROOM_SIZE,
    SPATIAL_GRID_DIMENSION,
    THREADS,
    MAX_ENTITIES,
    MAX_ENTITIES_TO_REPLACE,
    COLLISION_SEPARATION_CONSTANT,
    STORE_COLLISIONS
};
use std::{thread, ptr, collections::HashSet, sync::Mutex};
use wide::f32x4;

// DO NOT SET THESE
const ENCODING_BITS: usize = (usize::BITS - SPATIAL_GRID_DIMENSION.leading_zeros() - 1) as usize;
const SPATIAL_GRID_AREA: usize = SPATIAL_GRID_DIMENSION * SPATIAL_GRID_DIMENSION;
const ROOM_GRID_RATIO: f32 = (SPATIAL_GRID_DIMENSION as f32) / ROOM_SIZE;

#[derive(Clone)]
pub struct Entity {
    pub index: usize,
    grid_pos_x: usize,
    grid_pos_y: usize,
    grid_body: usize,
    pub x: f32,
    pub y: f32,
    pub mass: f32,
    pub inv_mass: f32,
    pub friction: f32,
    pub velocity_x: f32,
    pub velocity_y: f32,
    pub terminal_velocity_in_direction: f32,
    pub acceleration_x: f32,
    pub acceleration_y: f32,
    pub movement_acceleration: f32,
    pub movement_acceleration_45_deg: f32,
    pub radius: f32,
    pub body_type: u8,
    has_friction: bool,
    pub replace: bool
}

pub struct Room {
    pub spatial_grid: Vec<Mutex<HashSet<usize>>>,
    pub spatial_grid_ptr: *mut Mutex<HashSet<usize>>,
    pub entities: Vec<Entity>,
    pub entities_ptr: *mut Entity,
    applied_collisions: Vec<[f32; 4]>,
    applied_collisions_ptr: *mut [f32; 4],
    pub entity_chunks: [[usize; 2]; THREADS],
    pub collision_chunks: [[usize; 2]; THREADS],
    pub replacement_queue: Vec<usize>,
    pub replacement_queue_ptr: *mut usize,
    pub stored_collisions: Vec<Vec<[usize; 2]>>,
    pub stored_collisions_ptrs: [*mut Vec<[usize; 2]>; THREADS],
    pub tick: usize
}

#[inline(always)]
fn manage_collision_circle_circle(dsq: f32, entity: &Entity, collision_entity: &Entity, applied_collisions_ptrs_ptr: *mut [f32; 4], thread: usize, entities_size: usize, stored_collisions: &mut Vec<[usize; 2]>) {
    let distance = dsq.sqrt();
    if distance == 0.0 { return };
    let inv_distance = 1.0 / distance;
    let nx = (collision_entity.x - entity.x) * inv_distance;
    let ny = (collision_entity.y - entity.y) * inv_distance;
    let separation = (entity.radius + collision_entity.radius) - distance + COLLISION_SEPARATION_CONSTANT;
    if separation <= 0.0 { return };
    let inv_total_mass = 1.0 / (entity.mass + collision_entity.mass);
    let ratio_entity = collision_entity.mass * inv_total_mass;
    let ratio_collision = entity.mass * inv_total_mass;
    let v1n = entity.velocity_x * nx + entity.velocity_y * ny;
    let v2n = collision_entity.velocity_x * nx + collision_entity.velocity_y * ny;
    let v1n_new = ((entity.mass - collision_entity.mass) * v1n + 2.0 * collision_entity.mass * v2n) * inv_total_mass;
    let v2n_new = ((collision_entity.mass - entity.mass) * v2n + 2.0 * entity.mass * v1n) * inv_total_mass;
    let dv1 = v1n_new - v1n;
    let dv2 = v2n_new - v2n;
    let data_entity = [dv1 * nx, dv1 * ny, -nx * separation * ratio_entity, -ny * separation * ratio_entity];
    let data_collision_entity = [dv2 * nx, dv2 * ny, nx * separation * ratio_collision, ny * separation * ratio_collision];
    unsafe { 
        let entity_pos = applied_collisions_ptrs_ptr.add(thread * entities_size + entity.index);
        let collision_entity_pos = applied_collisions_ptrs_ptr.add(thread * entities_size + collision_entity.index);
        *entity_pos = (f32x4::from(*entity_pos) + f32x4::from(data_entity)).to_array();
        *collision_entity_pos = (f32x4::from(*collision_entity_pos) + f32x4::from(data_collision_entity)).to_array();
    };
    if STORE_COLLISIONS {
        stored_collisions.push([entity.index, collision_entity.index]);
    };
}

#[inline(always)]
fn manage_collision_circle_square(cx: f32, cy: f32, circle: &Entity, square: &Entity, applied_collisions_ptrs_ptr: *mut [f32; 4], thread: usize, entities_size: usize, stored_collisions: &mut Vec<[usize; 2]>) {
    let dx = circle.x - cx;
    let dy = circle.y - cy;
    let dist_sq = dx * dx + dy * dy;
    let (nx, ny, separation) = if dist_sq == 0.0 {
        let d_left = (circle.x - (square.x - square.radius)).abs();
        let d_right = (circle.x - (square.x + square.radius)).abs();
        let d_top = (circle.y - (square.y - square.radius)).abs();
        let d_bottom = (circle.y - (square.y + square.radius)).abs();
        let min_dist = d_left.min(d_right).min(d_top).min(d_bottom);
        if min_dist == d_left { (-1.0, 0.0, circle.radius + min_dist + COLLISION_SEPARATION_CONSTANT) }
        else if min_dist == d_right { (1.0, 0.0, circle.radius + min_dist + COLLISION_SEPARATION_CONSTANT) }
        else if min_dist == d_top { (0.0, -1.0, circle.radius + min_dist + COLLISION_SEPARATION_CONSTANT) }
        else { (0.0, 1.0, circle.radius + min_dist + COLLISION_SEPARATION_CONSTANT) }
    } else {
        let dist = dist_sq.sqrt();
        let inv_dist = 1.0 / dist;
        (dx * inv_dist, dy * inv_dist, circle.radius - dist + COLLISION_SEPARATION_CONSTANT)
    };
    let inv_total_mass = 1.0 / (circle.mass + square.mass);
    let ratio_circle = square.mass * inv_total_mass;
    let ratio_square = circle.mass * inv_total_mass;
    let dvx = circle.velocity_x - square.velocity_x;
    let dvy = circle.velocity_y - square.velocity_y;
    let vn = dvx * nx + dvy * ny;
    let (impulse_cx, impulse_cy, impulse_sx, impulse_sy) = if vn < 0.0 {
        let impulse = (2.0 * vn) * inv_total_mass;
        (-impulse * square.mass * nx, -impulse * square.mass * ny, impulse * circle.mass * nx, impulse * circle.mass * ny)
    } else { (0.0, 0.0, 0.0, 0.0) };
    let data_circle = [impulse_cx, impulse_cy, nx * separation * ratio_circle, ny * separation * ratio_circle];
    let data_square = [impulse_sx, impulse_sy, -nx * separation * ratio_square, -ny * separation * ratio_square];
    unsafe {
        let circle_pos = applied_collisions_ptrs_ptr.add(thread * entities_size + circle.index);
        let square_pos = applied_collisions_ptrs_ptr.add(thread * entities_size + square.index);
        *circle_pos = (f32x4::from(*circle_pos) + f32x4::from(data_circle)).to_array();
        *square_pos = (f32x4::from(*square_pos) + f32x4::from(data_square)).to_array();
    };
    if STORE_COLLISIONS {
        stored_collisions.push([circle.index, square.index]);
    };
}

#[inline(always)]
fn manage_collision_square_square(entity: &Entity, collision_entity: &Entity, applied_collisions_ptrs_ptr: *mut [f32; 4], thread: usize, entities_size: usize, stored_collisions: &mut Vec<[usize; 2]>) {
    let overlap_x = entity.radius + collision_entity.radius - f32::abs(collision_entity.x - entity.x);
    let overlap_y = entity.radius + collision_entity.radius - f32::abs(collision_entity.y - entity.y);
    if overlap_x <= 0.0 || overlap_y <= 0.0 { return };
    let (nx, ny, separation) = if overlap_x < overlap_y {
        (if collision_entity.x - entity.x > 0.0 { 1.0 } else { -1.0 }, 0.0, overlap_x + COLLISION_SEPARATION_CONSTANT)
    } else {
        (0.0, if collision_entity.y - entity.y > 0.0 { 1.0 } else { -1.0 }, overlap_y + COLLISION_SEPARATION_CONSTANT)
    };
    let inv_m1 = if entity.mass == f32::MAX { 0.0 } else { entity.inv_mass };
    let inv_m2 = if collision_entity.mass == f32::MAX { 0.0 } else { collision_entity.inv_mass };
    let inv_mass_sum = inv_m1 + inv_m2;
    let inv_sum_final = 1.0 / inv_mass_sum;
    let v1n = entity.velocity_x * nx + entity.velocity_y * ny;
    let v2n = collision_entity.velocity_x * nx + collision_entity.velocity_y * ny;
    let rel_vn = v1n - v2n;
    let j = -(2.0 * rel_vn) * inv_sum_final;
    let data_entity = [j * inv_m1 * nx, j * inv_m1 * ny, -nx * separation * (inv_m1 * inv_sum_final), -ny * separation * (inv_m1 * inv_sum_final)];
    let data_collision_entity = [-j * inv_m2 * nx, -j * inv_m2 * ny, nx * separation * (inv_m2 * inv_sum_final), ny * separation * (inv_m2 * inv_sum_final)];
    unsafe { 
        let entity_pos = applied_collisions_ptrs_ptr.add(thread * entities_size + entity.index);
        let collision_entity_pos = applied_collisions_ptrs_ptr.add(thread * entities_size + collision_entity.index);
        *entity_pos = (f32x4::from(*entity_pos) + f32x4::from(data_entity)).to_array();
        *collision_entity_pos = (f32x4::from(*collision_entity_pos) + f32x4::from(data_collision_entity)).to_array();
    };
    if STORE_COLLISIONS {
        stored_collisions.push([entity.index, collision_entity.index]);
    };
}

#[inline(always)]
fn update_grid_position(pos: usize, index: usize, spatial_grid: *mut Mutex<HashSet<usize>>) {
    unsafe { 
        let mut position = (*spatial_grid.add(pos)).lock().unwrap();
        position.insert(index); 
    };
}

#[inline(always)]
fn remove_grid_position(pos: usize, index: &usize, spatial_grid: *mut Mutex<HashSet<usize>>) {
    unsafe { 
        let mut position = (*spatial_grid.add(pos)).lock().unwrap();
        position.remove(index); 
    };
}

#[inline(always)]
fn update_entity_body(y_bound: usize, x_bound: usize, y_position_update: usize, x_position_update: usize, y_position_delete: usize, x_position_delete: usize, index: usize, spatial_grid: *mut Mutex<HashSet<usize>>) {
    for y_offset in y_position_update..usize::min(y_position_update + y_bound, SPATIAL_GRID_DIMENSION) {
        for x_offset in x_position_update..usize::min(x_position_update + x_bound, SPATIAL_GRID_DIMENSION) {
            let position = (y_offset << ENCODING_BITS) | x_offset;
            update_grid_position(position, index, spatial_grid);
        };
    };
    for y_offset in y_position_delete..usize::min(y_position_delete + y_bound, SPATIAL_GRID_DIMENSION) {
        for x_offset in x_position_delete..usize::min(x_position_delete + x_bound, SPATIAL_GRID_DIMENSION) {
            let dposition = (y_offset << ENCODING_BITS) | x_offset;
            remove_grid_position(dposition, &index, spatial_grid);
        };
    };
}

impl Room {
    #[inline(always)]
    pub fn init() -> Room {
        let mut entities = Vec::with_capacity(MAX_ENTITIES);
        let entities_ptr = entities.as_mut_ptr();
        let mut replacement_queue = Vec::with_capacity(MAX_ENTITIES_TO_REPLACE);
        let replacement_queue_ptr = replacement_queue.as_mut_ptr();
        let mut spatial_grid = Vec::with_capacity(SPATIAL_GRID_AREA);
        for _pos in 0..SPATIAL_GRID_AREA {
            spatial_grid.push(Mutex::new(HashSet::new()));
        };
        let spatial_grid_ptr = spatial_grid.as_mut_ptr();
        let mut collision_chunks = [[0, 0]; THREADS];
        let chunk_size = SPATIAL_GRID_AREA / THREADS;
        let chunk_rem = SPATIAL_GRID_AREA % THREADS;
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
        let mut applied_collisions = Vec::with_capacity(THREADS * MAX_ENTITIES);
        let applied_collisions_ptr = applied_collisions.as_mut_ptr();
        let mut stored_collisions = Vec::with_capacity(THREADS);
        let stored_collisions_ptr = stored_collisions.as_mut_ptr();
        let mut stored_collisions_ptrs = [ptr::null_mut(); THREADS];
        unsafe { stored_collisions.set_len(THREADS) };
        for bucket in 0..THREADS {
            stored_collisions[bucket] = Vec::new();
            unsafe { stored_collisions_ptrs[bucket] = stored_collisions_ptr.add(bucket) };
        };
        Room {
            spatial_grid,
            spatial_grid_ptr,
            entities,
            entities_ptr,
            applied_collisions,
            applied_collisions_ptr,
            entity_chunks: [[0, 0]; THREADS],
            collision_chunks,
            replacement_queue,
            replacement_queue_ptr,
            stored_collisions,
            stored_collisions_ptrs,
            tick: 0
        }
    }

    #[inline(always)]
    fn update_chunks(&mut self) {
        let len = self.entities.len();
        let chunk_size = len / THREADS;
        let chunk_rem = len % THREADS;
        let mut chunk_pos = 0;
        for chunk in 0..THREADS {
            let next_chunk = if chunk < chunk_rem {
                chunk_pos + chunk_size + 1
            } else {
                chunk_pos + chunk_size
            };
            self.entity_chunks[chunk] = [chunk_pos, next_chunk];
            chunk_pos = next_chunk;
        };
    }

    #[inline(always)]
    pub fn create_entity(&mut self, x: f32, y: f32, mass: f32, friction: f32, velocity_x: f32, velocity_y: f32, terminal_velocity_in_direction: f32, movement_acceleration: f32, radius: f32, body_type: u8, has_friction: bool) -> usize {
        let grid_body = (radius * 2.0 * ROOM_GRID_RATIO).ceil() as usize + 1;
        let grid_pos_x = ((x - radius) * ROOM_GRID_RATIO) as usize;
        let grid_pos_y = ((y - radius) * ROOM_GRID_RATIO) as usize;
        let index: usize;
        if self.replacement_queue.len() == 0 {
            unsafe {
                index = self.entities.len();
                if self.entities.len() == self.entities.capacity() {
                    let old_len = self.entities.len();
                    let new_capacity = self.entities.capacity() * 2;
                    let mut new_entities = Vec::with_capacity(new_capacity);
                    ptr::copy_nonoverlapping(self.entities_ptr, new_entities.as_mut_ptr(), old_len);
                    new_entities.set_len(old_len);
                    self.entities = new_entities;
                    self.entities_ptr = self.entities.as_mut_ptr();
                };
                *self.entities_ptr.add(self.entities.len()) = Entity {
                    index,
                    grid_pos_x,
                    grid_pos_y,
                    grid_body,
                    x,
                    y,
                    mass,
                    inv_mass: 1.0 / mass,
                    friction,
                    velocity_x,
                    velocity_y,
                    terminal_velocity_in_direction,
                    acceleration_x: 0.0,
                    acceleration_y: 0.0,
                    movement_acceleration,
                    movement_acceleration_45_deg: movement_acceleration * 0.70710675,
                    radius,
                    body_type,
                    has_friction,
                    replace: false
                };
                self.entities.set_len(self.entities.len() + 1);
                self.update_chunks();
                if self.entities.len() * THREADS > self.applied_collisions.capacity() {
                    self.applied_collisions = Vec::with_capacity(self.applied_collisions.capacity() * 2);
                    self.applied_collisions_ptr = self.applied_collisions.as_mut_ptr();
                };
            };
        } else {
            unsafe { 
                let len = self.replacement_queue.len() - 1;
                index = ptr::read(self.replacement_queue_ptr.add(len));
                self.replacement_queue.set_len(len);
                *self.entities_ptr.add(index) = Entity {
                    index,
                    grid_pos_x,
                    grid_pos_y,
                    grid_body,
                    x,
                    y,
                    mass,
                    inv_mass: 1.0 / mass,
                    friction,
                    velocity_x,
                    velocity_y,
                    terminal_velocity_in_direction,
                    acceleration_x: 0.0,
                    acceleration_y: 0.0,
                    movement_acceleration,
                    movement_acceleration_45_deg: movement_acceleration * 0.70710675,
                    radius,
                    body_type,
                    has_friction,
                    replace: false
                }; 
            };
        };
        for y_pos in grid_pos_y..usize::min(grid_pos_y + grid_body, SPATIAL_GRID_DIMENSION) {
            for x_pos in grid_pos_x..usize::min(grid_pos_x + grid_body, SPATIAL_GRID_DIMENSION) {
                update_grid_position((y_pos << ENCODING_BITS) | x_pos, index, self.spatial_grid_ptr);
            };
        };
        index
    }

    #[inline(always)]
    pub fn remove_entity(&mut self, index: usize) {
        unsafe { 
            let entity = self.entities.get_unchecked_mut(index);
            if entity.replace { return };
            entity.replace = true;
            if self.replacement_queue.len() == self.replacement_queue.capacity() {
                let old_len = self.replacement_queue.len();
                let new_capacity = self.replacement_queue.capacity() * 2;
                let mut new_queue = Vec::with_capacity(new_capacity);
                ptr::copy_nonoverlapping(self.replacement_queue_ptr, new_queue.as_mut_ptr(), old_len);
                new_queue.set_len(old_len);
                let new_queue_ptr = new_queue.as_mut_ptr();
                self.replacement_queue = new_queue;
                self.replacement_queue_ptr = new_queue_ptr;
            };
            *self.replacement_queue_ptr.add(self.replacement_queue.len()) = entity.index;
            self.replacement_queue.set_len(self.replacement_queue.len() + 1);
            for y_pos in entity.grid_pos_y..usize::min(entity.grid_pos_y + entity.grid_body, SPATIAL_GRID_DIMENSION) {
                for x_pos in entity.grid_pos_x..usize::min(entity.grid_pos_x + entity.grid_body, SPATIAL_GRID_DIMENSION) {
                    remove_grid_position((y_pos << ENCODING_BITS) | x_pos, &entity.index, self.spatial_grid_ptr);
                };
            };
        };
    }

    #[inline(always)]
    pub fn create_entity_movement_from_angle(&mut self, index: usize, angle: f32) {
        let entity = unsafe { self.entities.get_unchecked_mut(index) };
        entity.acceleration_x = f32::cos(angle) * entity.movement_acceleration;
        entity.acceleration_y = f32::sin(angle) * entity.movement_acceleration;
    }

    #[inline(always)]
    pub fn create_entity_movement_from_cardinal_direction(&mut self, index: usize, direction: u8) {
        let entity = unsafe { self.entities.get_unchecked_mut(index) };
        match direction {
            0 => {
                entity.acceleration_x = entity.movement_acceleration;
            }
            1 => {
                entity.acceleration_x = -entity.movement_acceleration;
            }
            2 => {
                entity.acceleration_y = entity.movement_acceleration;
            }
            3 => {
                entity.acceleration_y = -entity.movement_acceleration;
            }
            4 => {
                entity.acceleration_x = entity.movement_acceleration_45_deg;
                entity.acceleration_y = entity.movement_acceleration_45_deg;
            }
            5 => {
                entity.acceleration_x = entity.movement_acceleration_45_deg;
                entity.acceleration_y = -entity.movement_acceleration_45_deg;
            }
            6 => {
                entity.acceleration_x = -entity.movement_acceleration_45_deg;
                entity.acceleration_y = entity.movement_acceleration_45_deg;
            }
            7 => {
                entity.acceleration_x = -entity.movement_acceleration_45_deg;
                entity.acceleration_y = -entity.movement_acceleration_45_deg;
            }
            _ => {}
        }
    }

    #[inline(always)]
    pub fn stop_entity_movement(&mut self, index: usize) {
        let entity = unsafe { self.entities.get_unchecked_mut(index) };
        entity.acceleration_x = 0.0;
        entity.acceleration_y = 0.0;
    }

    #[inline(always)]
    pub fn stop_entity_movement_x(&mut self, index: usize) {
        let entity = unsafe { self.entities.get_unchecked_mut(index) };
        entity.acceleration_x = 0.0;
    }

    #[inline(always)]
    pub fn stop_entity_movement_y(&mut self, index: usize) {
        let entity = unsafe { self.entities.get_unchecked_mut(index) };
        entity.acceleration_y = 0.0;
    }

    #[inline(always)]
    pub fn update(&mut self) {
        let entities_ptr = self.entities_ptr as usize;
        let spatial_grid_ptr = self.spatial_grid_ptr as usize;
        let mut spatial_grid_update_threads = Vec::with_capacity(THREADS);
        for thread in 0..THREADS {
            let chunk = self.entity_chunks[thread];
            let worker = thread::spawn(move || {
                let entities = entities_ptr as *mut Entity;
                let spatial_grid = spatial_grid_ptr as *mut Mutex<HashSet<usize>>;
                for index in chunk[0]..chunk[1] {
                    let entity = unsafe { &mut *entities.add(index) };
                    if !entity.replace {
                        entity.velocity_x += entity.acceleration_x;
                        entity.velocity_y += entity.acceleration_y;
                        entity.x += entity.velocity_x;
                        entity.y += entity.velocity_y;
                        if entity.x < entity.radius { 
                            entity.x = entity.radius; 
                            entity.velocity_x = -entity.velocity_x; 
                        } else { 
                            if entity.x > ROOM_SIZE - entity.radius { 
                                entity.x = ROOM_SIZE - entity.radius; 
                                entity.velocity_x = -entity.velocity_x; 
                            }; 
                        };
                        if entity.y < entity.radius { 
                            entity.y = entity.radius; 
                            entity.velocity_y = -entity.velocity_y; 
                        } else { 
                            if entity.y > ROOM_SIZE - entity.radius { 
                                entity.y = ROOM_SIZE - entity.radius; 
                                entity.velocity_y = -entity.velocity_y; 
                            }; 
                        };
                        if entity.has_friction {
                            entity.velocity_x *= entity.friction;
                            entity.velocity_y *= entity.friction;
                        };
                        let speed_sq = entity.velocity_x * entity.velocity_x + entity.velocity_y * entity.velocity_y;
                        if speed_sq > entity.terminal_velocity_in_direction * entity.terminal_velocity_in_direction { 
                            let speed_scalar = entity.terminal_velocity_in_direction / speed_sq.sqrt(); 
                            entity.velocity_x *= speed_scalar; 
                            entity.velocity_y *= speed_scalar; 
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
                                update_entity_body(entity.grid_body, entity.grid_body, spatial_grid_pos_y, spatial_grid_pos_x, entity.grid_pos_y, entity.grid_pos_x, entity.index, spatial_grid);
                            } else {
                                if shift_x.0 == 0 {
                                    if shift_y.1 {
                                        update_entity_body(shift_y.0, entity.grid_body, entity.grid_pos_y + entity.grid_body, spatial_grid_pos_x, entity.grid_pos_y, spatial_grid_pos_x, entity.index, spatial_grid);
                                    } else {
                                        update_entity_body(shift_y.0, entity.grid_body, spatial_grid_pos_y, spatial_grid_pos_x, spatial_grid_pos_y + entity.grid_body, spatial_grid_pos_x, entity.index, spatial_grid);
                                    };
                                } else if shift_y.0 == 0 {
                                    if shift_x.1 {
                                        update_entity_body(entity.grid_body, shift_x.0, spatial_grid_pos_y, entity.grid_pos_x + entity.grid_body, spatial_grid_pos_y, entity.grid_pos_x, entity.index, spatial_grid);
                                    } else {
                                        update_entity_body(entity.grid_body, shift_x.0, spatial_grid_pos_y, spatial_grid_pos_x, spatial_grid_pos_y, spatial_grid_pos_x + entity.grid_body, entity.index, spatial_grid);
                                    };
                                } else {
                                    if shift_y.1 {
                                        update_entity_body(shift_y.0, entity.grid_body, spatial_grid_pos_y + entity.grid_body - shift_y.0, spatial_grid_pos_x, entity.grid_pos_y, entity.grid_pos_x, entity.index, spatial_grid);
                                    } else {
                                        update_entity_body(shift_y.0, entity.grid_body, spatial_grid_pos_y, spatial_grid_pos_x, entity.grid_pos_y + entity.grid_body - shift_y.0, entity.grid_pos_x, entity.index, spatial_grid);
                                    };
                                    if shift_x.1 {
                                        update_entity_body(entity.grid_body - shift_y.0, shift_x.0, spatial_grid_pos_y, spatial_grid_pos_x + entity.grid_body - shift_x.0, entity.grid_pos_y, entity.grid_pos_x, entity.index, spatial_grid);
                                    } else {
                                        update_entity_body(entity.grid_body - shift_y.0, shift_x.0, spatial_grid_pos_y, spatial_grid_pos_x, entity.grid_pos_y, entity.grid_pos_x + entity.grid_body - shift_x.0, entity.index, spatial_grid);
                                    };
                                };
                            };
                            entity.grid_pos_x = spatial_grid_pos_x;
                            entity.grid_pos_y = spatial_grid_pos_y; 
                        };
                    };
                };
            });
            spatial_grid_update_threads.push(worker);
        };
        for thread in spatial_grid_update_threads {
            thread.join().unwrap();
        };
        let entities_size = self.entities.len();
        let mut collision_storer_initializer_threads = Vec::with_capacity(THREADS);
        unsafe { self.applied_collisions.set_len(THREADS * entities_size) };
        let applied_collisions_ptr = self.applied_collisions_ptr as usize;
        for thread in 0..THREADS {
            let worker = thread::spawn(move || {
                let applied_collisions = applied_collisions_ptr as *mut [f32; 4];
                unsafe { 
                    let bucket_start = applied_collisions.add(thread * entities_size);
                    ptr::write_bytes(bucket_start, 0, entities_size);
                };
            });
            collision_storer_initializer_threads.push(worker);
        };
        for worker in collision_storer_initializer_threads {
            worker.join().unwrap();
        };
        let mut collision_threads = Vec::with_capacity(THREADS);
        let stored_collisions_ptr = self.stored_collisions_ptrs.as_mut_ptr() as usize;
        for thread in 0..THREADS {
            let chunk = self.collision_chunks[thread];
            let worker = thread::spawn(move || {
                let applied_collisions = applied_collisions_ptr as *mut [f32; 4];
                let entities = entities_ptr as *mut Entity;
                let spatial_grid = spatial_grid_ptr as *mut Mutex<HashSet<usize>>;
                let stored_collisions = unsafe { &mut **(stored_collisions_ptr as *mut *mut Vec<[usize; 2]>).add(thread) };
                stored_collisions.clear();
                for grid_index in chunk[0]..chunk[1] {
                    let spatial_grid_position = unsafe { (*spatial_grid.add(grid_index)).lock().unwrap() };
                    if spatial_grid_position.len() == 0 || spatial_grid_position.len() == 1 { continue };
                    let entities_vec: Vec<usize> = spatial_grid_position.iter().copied().collect();
                    let entities_vec_ptr = entities_vec.as_ptr();
                    for entity_index in 0..entities_vec.len() {
                        for collision_index in entity_index + 1..entities_vec.len() {
                            let entity = unsafe { &*entities.add(*entities_vec_ptr.add(entity_index)) };
                            let collision_entity = unsafe { &*entities.add(*entities_vec_ptr.add(collision_index)) };
                            let min_shared_x = usize::max(entity.grid_pos_x, collision_entity.grid_pos_x);
                            let min_shared_y = usize::max(entity.grid_pos_y, collision_entity.grid_pos_y);
                            let min_shared_coordinate = (min_shared_y << ENCODING_BITS) | min_shared_x;
                            if min_shared_coordinate == grid_index {
                                if entity.body_type + collision_entity.body_type == 2 {
                                    let dsq = (collision_entity.x - entity.x) * (collision_entity.x - entity.x) + (collision_entity.y - entity.y) * (collision_entity.y - entity.y);
                                    if dsq < (entity.radius + collision_entity.radius) * (entity.radius + collision_entity.radius) {
                                        manage_collision_circle_circle(dsq, entity, collision_entity, applied_collisions, thread, entities_size, stored_collisions);
                                    };
                                } else if entity.body_type + collision_entity.body_type == 1 {
                                    if entity.body_type == 0 {
                                        let cx = collision_entity.x.clamp(entity.x - entity.radius, entity.x + entity.radius);
                                        let cy = collision_entity.y.clamp(entity.y - entity.radius, entity.y + entity.radius);
                                        if (cx - collision_entity.x) * (cx - collision_entity.x) + (cy - collision_entity.y) * (cy - collision_entity.y) < collision_entity.radius * collision_entity.radius {
                                            manage_collision_circle_square(cx, cy, collision_entity, entity, applied_collisions, thread, entities_size, stored_collisions);
                                        };
                                    } else {
                                        let cx = entity.x.clamp(collision_entity.x - collision_entity.radius, collision_entity.x + collision_entity.radius);
                                        let cy = entity.y.clamp(collision_entity.y - collision_entity.radius, collision_entity.y + collision_entity.radius);
                                        if (cx - entity.x) * (cx - entity.x) + (cy - entity.y) * (cy - entity.y) < entity.radius * entity.radius {
                                            manage_collision_circle_square(cx, cy, entity, collision_entity, applied_collisions, thread, entities_size, stored_collisions);
                                        };
                                    };
                                } else {
                                    if f32::abs(collision_entity.x - entity.x) < entity.radius + collision_entity.radius && f32::abs(collision_entity.y - entity.y) < entity.radius + collision_entity.radius {
                                        manage_collision_square_square(entity, collision_entity, applied_collisions, thread, entities_size, stored_collisions);
                                    };
                                };
                            };
                        };
                    };
                };
            });
            collision_threads.push(worker);
        };
        for thread in collision_threads {
            thread.join().unwrap();
        };
        let mut collision_result_handler_threads = Vec::with_capacity(THREADS);
        for thread in 0..THREADS {
            let chunk = self.entity_chunks[thread];
            let worker = thread::spawn(move || {
                let applied_collisions = applied_collisions_ptr as *mut [f32; 4];
                let entities = entities_ptr as *mut Entity;
                unsafe { 
                    for index in chunk[0]..chunk[1] {
                        let entity = &mut *entities.add(index);
                        if entity.replace { continue };
                        let mut total_simd = f32x4::ZERO;
                        for bucket in 0..THREADS {
                            let bucket_pos = applied_collisions.add(bucket * entities_size + index);
                            total_simd += f32x4::from(*bucket_pos);
                        };
                        let total = total_simd.to_array();
                        entity.velocity_x += total[0];
                        entity.velocity_y += total[1];
                        entity.x += total[2];
                        entity.y += total[3];
                    };
                };
            });
            collision_result_handler_threads.push(worker);
        };
        for worker in collision_result_handler_threads {
            worker.join().unwrap();
        };
        self.tick += 1;
    }
}
