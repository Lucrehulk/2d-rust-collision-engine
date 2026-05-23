use super::config::config_data::{
    COLLISION_SEPARATION_CONSTANT, MAX_ENTITIES, MAX_ENTITIES_TO_REPLACE, CELL_SIZE,
    STORE_COLLISIONS, THREADS,
};
use rayon;
use std::{
    collections::HashSet,
    ptr,
    sync::atomic::{AtomicBool, Ordering},
};
use wide::f32x4;

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
    pub radius: f32,
    pub body_type: u8,
    pub replace: bool,
    pub collision_type: u8,
}

pub struct Room {
    pub spatial_grid: Vec<HashSet<usize>>,
    pub spatial_grid_locks: Vec<AtomicBool>,
    pub spatial_grid_dimension: usize,
    pub encoding_bits: usize,
    pub size: f32,
    pub grid_ratio: f32,
    pub entities: Vec<Entity>,
    applied_collisions: Vec<[f32; 4]>,
    pub entity_chunks: [[usize; 2]; THREADS],
    pub replacement_queue: Vec<usize>,
    pub stored_collisions: Vec<Vec<[usize; 2]>>,
    pub tick: usize,
}

#[inline(always)]
fn manage_collision_circle_circle(
    dsq: f32,
    entity: &Entity,
    collision_entity: &Entity,
    applied_collisions_ptrs_ptr: *mut [f32; 4],
    thread: usize,
    entities_size: usize,
    stored_collisions: &mut Vec<[usize; 2]>,
) {
    if (entity.collision_type == collision_entity.collision_type && entity.collision_type != 0)
        || (entity.collision_type < 3
            && collision_entity.collision_type < 3
            && entity.collision_type + collision_entity.collision_type > 1)
    {
        let distance = dsq.sqrt();
        if distance == 0.0 {
            return;
        };
        let inv_distance = 1.0 / distance;
        let nx = (collision_entity.x - entity.x) * inv_distance;
        let ny = (collision_entity.y - entity.y) * inv_distance;
        let separation =
            (entity.radius + collision_entity.radius) - distance + COLLISION_SEPARATION_CONSTANT;
        if separation <= 0.0 {
            return;
        };
        let inv_total_mass = 1.0 / (entity.mass + collision_entity.mass);
        let ratio_entity = collision_entity.mass * inv_total_mass;
        let ratio_collision = entity.mass * inv_total_mass;
        let v1n = entity.velocity_x * nx + entity.velocity_y * ny;
        let v2n = collision_entity.velocity_x * nx + collision_entity.velocity_y * ny;
        let (dv1, dv2) = if v1n - v2n > 0.0 {
            let v1n_new = ((entity.mass - collision_entity.mass) * v1n
                + 2.0 * collision_entity.mass * v2n)
                * inv_total_mass;
            let v2n_new = ((collision_entity.mass - entity.mass) * v2n + 2.0 * entity.mass * v1n)
                * inv_total_mass;
            (v1n_new - v1n, v2n_new - v2n)
        } else {
            (0.0, 0.0)
        };
        let data_entity = [
            dv1 * nx,
            dv1 * ny,
            -nx * separation * ratio_entity,
            -ny * separation * ratio_entity,
        ];
        let data_collision_entity = [
            dv2 * nx,
            dv2 * ny,
            nx * separation * ratio_collision,
            ny * separation * ratio_collision,
        ];
        unsafe {
            let entity_pos = applied_collisions_ptrs_ptr.add(thread * entities_size + entity.index);
            let collision_entity_pos =
                applied_collisions_ptrs_ptr.add(thread * entities_size + collision_entity.index);
            *entity_pos = (f32x4::from(*entity_pos) + f32x4::from(data_entity)).to_array();
            *collision_entity_pos = (f32x4::from(*collision_entity_pos)
                + f32x4::from(data_collision_entity))
            .to_array();
        };
    };
    if STORE_COLLISIONS {
        stored_collisions.push([entity.index, collision_entity.index]);
    };
}

#[inline(always)]
fn manage_collision_circle_square(
    cx: f32,
    cy: f32,
    circle: &Entity,
    square: &Entity,
    applied_collisions_ptrs_ptr: *mut [f32; 4],
    thread: usize,
    entities_size: usize,
    stored_collisions: &mut Vec<[usize; 2]>,
) {
    if (circle.collision_type == square.collision_type && circle.collision_type != 0)
        || (circle.collision_type < 3
            && square.collision_type < 3
            && circle.collision_type + square.collision_type > 1)
    {
        let dx = circle.x - cx;
        let dy = circle.y - cy;
        let dist_sq = dx * dx + dy * dy;
        let (nx, ny, separation) = if dist_sq == 0.0 {
            let d_left = (circle.x - (square.x - square.radius)).abs();
            let d_right = (circle.x - (square.x + square.radius)).abs();
            let d_top = (circle.y - (square.y - square.radius)).abs();
            let d_bottom = (circle.y - (square.y + square.radius)).abs();
            let min_dist = d_left.min(d_right).min(d_top).min(d_bottom);
            if min_dist == d_left {
                (
                    -1.0,
                    0.0,
                    circle.radius + min_dist + COLLISION_SEPARATION_CONSTANT,
                )
            } else if min_dist == d_right {
                (
                    1.0,
                    0.0,
                    circle.radius + min_dist + COLLISION_SEPARATION_CONSTANT,
                )
            } else if min_dist == d_top {
                (
                    0.0,
                    -1.0,
                    circle.radius + min_dist + COLLISION_SEPARATION_CONSTANT,
                )
            } else {
                (
                    0.0,
                    1.0,
                    circle.radius + min_dist + COLLISION_SEPARATION_CONSTANT,
                )
            }
        } else {
            let dist = dist_sq.sqrt();
            let inv_dist = 1.0 / dist;
            (
                dx * inv_dist,
                dy * inv_dist,
                circle.radius - dist + COLLISION_SEPARATION_CONSTANT,
            )
        };
        let inv_total_mass = 1.0 / (circle.mass + square.mass);
        let ratio_circle = square.mass * inv_total_mass;
        let ratio_square = circle.mass * inv_total_mass;
        let dvx = circle.velocity_x - square.velocity_x;
        let dvy = circle.velocity_y - square.velocity_y;
        let vn = dvx * nx + dvy * ny;
        let (impulse_cx, impulse_cy, impulse_sx, impulse_sy) = if vn < 0.0 {
            let impulse = (2.0 * vn) * inv_total_mass;
            (
                -impulse * square.mass * nx,
                -impulse * square.mass * ny,
                impulse * circle.mass * nx,
                impulse * circle.mass * ny,
            )
        } else {
            (0.0, 0.0, 0.0, 0.0)
        };
        let data_circle = [
            impulse_cx,
            impulse_cy,
            nx * separation * ratio_circle,
            ny * separation * ratio_circle,
        ];
        let data_square = [
            impulse_sx,
            impulse_sy,
            -nx * separation * ratio_square,
            -ny * separation * ratio_square,
        ];
        unsafe {
            let circle_pos = applied_collisions_ptrs_ptr.add(thread * entities_size + circle.index);
            let square_pos = applied_collisions_ptrs_ptr.add(thread * entities_size + square.index);
            *circle_pos = (f32x4::from(*circle_pos) + f32x4::from(data_circle)).to_array();
            *square_pos = (f32x4::from(*square_pos) + f32x4::from(data_square)).to_array();
        };
    };
    if STORE_COLLISIONS {
        stored_collisions.push([circle.index, square.index]);
    };
}

#[inline(always)]
fn manage_collision_square_square(
    entity: &Entity,
    collision_entity: &Entity,
    applied_collisions_ptrs_ptr: *mut [f32; 4],
    thread: usize,
    entities_size: usize,
    stored_collisions: &mut Vec<[usize; 2]>,
) {
    if (entity.collision_type == collision_entity.collision_type && entity.collision_type != 0)
        || (entity.collision_type < 3
            && collision_entity.collision_type < 3
            && entity.collision_type + collision_entity.collision_type > 1)
    {
        let overlap_x =
            entity.radius + collision_entity.radius - f32::abs(collision_entity.x - entity.x);
        let overlap_y =
            entity.radius + collision_entity.radius - f32::abs(collision_entity.y - entity.y);
        if overlap_x <= 0.0 || overlap_y <= 0.0 {
            return;
        };
        let (nx, ny, separation) = if overlap_x < overlap_y {
            (
                if collision_entity.x - entity.x > 0.0 {
                    1.0
                } else {
                    -1.0
                },
                0.0,
                overlap_x + COLLISION_SEPARATION_CONSTANT,
            )
        } else {
            (
                0.0,
                if collision_entity.y - entity.y > 0.0 {
                    1.0
                } else {
                    -1.0
                },
                overlap_y + COLLISION_SEPARATION_CONSTANT,
            )
        };
        let inv_m1 = if entity.mass == f32::MAX {
            0.0
        } else {
            entity.inv_mass
        };
        let inv_m2 = if collision_entity.mass == f32::MAX {
            0.0
        } else {
            collision_entity.inv_mass
        };
        let inv_mass_sum = inv_m1 + inv_m2;
        let inv_sum_final = 1.0 / inv_mass_sum;
        let v1n = entity.velocity_x * nx + entity.velocity_y * ny;
        let v2n = collision_entity.velocity_x * nx + collision_entity.velocity_y * ny;
        let rel_vn = v1n - v2n;
        let j = if rel_vn > 0.0 {
            -(2.0 * rel_vn) * inv_sum_final
        } else {
            0.0
        };
        let data_entity = [
            j * inv_m1 * nx,
            j * inv_m1 * ny,
            -nx * separation * (inv_m1 * inv_sum_final),
            -ny * separation * (inv_m1 * inv_sum_final),
        ];
        let data_collision_entity = [
            -j * inv_m2 * nx,
            -j * inv_m2 * ny,
            nx * separation * (inv_m2 * inv_sum_final),
            ny * separation * (inv_m2 * inv_sum_final),
        ];
        unsafe {
            let entity_pos = applied_collisions_ptrs_ptr.add(thread * entities_size + entity.index);
            let collision_entity_pos =
                applied_collisions_ptrs_ptr.add(thread * entities_size + collision_entity.index);
            *entity_pos = (f32x4::from(*entity_pos) + f32x4::from(data_entity)).to_array();
            *collision_entity_pos = (f32x4::from(*collision_entity_pos)
                + f32x4::from(data_collision_entity))
            .to_array();
        };
    };
    if STORE_COLLISIONS {
        stored_collisions.push([entity.index, collision_entity.index]);
    };
}

#[inline(always)]
fn update_grid_position(
    pos: usize,
    index: usize,
    spatial_grid: *mut HashSet<usize>,
    spatial_grid_locks: *mut AtomicBool,
) {
    unsafe {
        while (*spatial_grid_locks.add(pos))
            .compare_exchange_weak(false, true, Ordering::Acquire, Ordering::Relaxed)
            .is_err()
        {
            std::hint::spin_loop();
        }
        (*spatial_grid.add(pos)).insert(index);
        (*spatial_grid_locks.add(pos)).store(false, Ordering::Release);
    };
}

#[inline(always)]
fn remove_grid_position(
    pos: usize,
    index: &usize,
    spatial_grid: *mut HashSet<usize>,
    spatial_grid_locks: *mut AtomicBool,
) {
    unsafe {
        while (*spatial_grid_locks.add(pos))
            .compare_exchange_weak(false, true, Ordering::Acquire, Ordering::Relaxed)
            .is_err()
        {
            std::hint::spin_loop();
        }
        (*spatial_grid.add(pos)).remove(index);
        (*spatial_grid_locks.add(pos)).store(false, Ordering::Release);
    };
}

#[inline(always)]
fn update_entity_body(
    y_bound: usize,
    x_bound: usize,
    y_position_update: usize,
    x_position_update: usize,
    y_position_delete: usize,
    x_position_delete: usize,
    index: usize,
    spatial_grid: *mut HashSet<usize>,
    spatial_grid_locks: *mut AtomicBool,
    spatial_grid_dimension: usize,
    encoding_bits: usize
) {
    for y_offset in
        y_position_update..usize::min(y_position_update + y_bound, spatial_grid_dimension)
    {
        for x_offset in
            x_position_update..usize::min(x_position_update + x_bound, spatial_grid_dimension)
        {
            let position = (y_offset << encoding_bits) | x_offset;
            update_grid_position(position, index, spatial_grid, spatial_grid_locks);
        }
    }
    for y_offset in
        y_position_delete..usize::min(y_position_delete + y_bound, spatial_grid_dimension)
    {
        for x_offset in
            x_position_delete..usize::min(x_position_delete + x_bound, spatial_grid_dimension)
        {
            let dposition = (y_offset << encoding_bits) | x_offset;
            remove_grid_position(dposition, &index, spatial_grid, spatial_grid_locks);
        }
    }
}

impl Room {
    #[inline(always)]
    pub fn init(size: f32) -> Room {
        let spatial_grid_dimension = (size / CELL_SIZE) as usize;
        let encoding_bits = (usize::BITS - spatial_grid_dimension.leading_zeros() - 1) as usize;
        let grid_ratio = (spatial_grid_dimension as f32) / size;
        let spatial_grid_area = spatial_grid_dimension * spatial_grid_dimension;

        let mut spatial_grid = Vec::with_capacity(spatial_grid_area);
        let mut spatial_grid_locks = Vec::with_capacity(spatial_grid_area);
        for _ in 0..spatial_grid_area {
            spatial_grid.push(HashSet::new());
            spatial_grid_locks.push(AtomicBool::new(false));
        }

        let mut stored_collisions = Vec::with_capacity(THREADS);
        for _ in 0..THREADS {
            stored_collisions.push(Vec::new());
        }

        Room {
            spatial_grid,
            spatial_grid_locks,
            spatial_grid_dimension,
            encoding_bits,
            size,
            grid_ratio,
            entities: Vec::with_capacity(MAX_ENTITIES),
            applied_collisions: Vec::with_capacity(THREADS * MAX_ENTITIES),
            entity_chunks: [[0, 0]; THREADS],
            replacement_queue: Vec::with_capacity(MAX_ENTITIES_TO_REPLACE),
            stored_collisions,
            tick: 0,
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
        }
    }

    #[inline(always)]
    pub fn create_entity(
        &mut self,
        mut x: f32,
        mut y: f32,
        mass: f32,
        friction: f32,
        velocity_x: f32,
        velocity_y: f32,
        terminal_velocity_in_direction: f32,
        movement_acceleration: f32,
        radius: f32,
        body_type: u8,
        collision_type: u8,
    ) -> usize {
        if x - radius < 0.0 {
            x = radius;
        } else if x + radius > self.size {
            x = self.size - radius;
        };
        if y - radius < 0.0 {
            y = radius;
        } else if y + radius > self.size {
            y = self.size - radius;
        };

        let grid_body = (radius * 2.0 * self.grid_ratio).ceil() as usize + 1;
        let grid_pos_x = ((x - radius) * self.grid_ratio) as usize;
        let grid_pos_y = ((y - radius) * self.grid_ratio) as usize;

        let index: usize;

        if let Some(reused_index) = self.replacement_queue.pop() {
            index = reused_index;
            self.entities[index] = Entity {
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
                radius,
                body_type,
                replace: false,
                collision_type,
            };
        } else {
            index = self.entities.len();
            self.entities.push(Entity {
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
                radius,
                body_type,
                replace: false,
                collision_type,
            });
            self.update_chunks();
        };

        let spatial_grid_ptr = self.spatial_grid.as_mut_ptr();
        let spatial_grid_locks_ptr = self.spatial_grid_locks.as_mut_ptr();

        for y_pos in grid_pos_y..usize::min(grid_pos_y + grid_body, self.spatial_grid_dimension) {
            for x_pos in grid_pos_x..usize::min(grid_pos_x + grid_body, self.spatial_grid_dimension) {
                update_grid_position(
                    (y_pos << self.encoding_bits) | x_pos,
                    index,
                    spatial_grid_ptr,
                    spatial_grid_locks_ptr,
                );
            }
        }
        index
    }

    #[inline(always)]
    pub fn resize_entity(&mut self, index: usize, size_change: f32, delta_size: bool) -> f32 {
        let spatial_grid_ptr = self.spatial_grid.as_mut_ptr();
        let spatial_grid_locks_ptr = self.spatial_grid_locks.as_mut_ptr();
        let entity = &mut self.entities[index];

        for y_pos in entity.grid_pos_y
            ..usize::min(entity.grid_pos_y + entity.grid_body, self.spatial_grid_dimension)
        {
            for x_pos in entity.grid_pos_x
                ..usize::min(entity.grid_pos_x + entity.grid_body, self.spatial_grid_dimension)
            {
                remove_grid_position(
                    (y_pos << self.encoding_bits) | x_pos,
                    &index,
                    spatial_grid_ptr,
                    spatial_grid_locks_ptr,
                );
            }
        }

        if delta_size {
            entity.radius += size_change;
        } else {
            entity.radius = size_change;
        };

        if entity.x - entity.radius < 0.0 {
            entity.x = entity.radius;
        } else if entity.x + entity.radius > self.size {
            entity.x = self.size - entity.radius;
        };
        if entity.y - entity.radius < 0.0 {
            entity.y = entity.radius;
        } else if entity.y + entity.radius > self.size {
            entity.y = self.size - entity.radius;
        };

        entity.grid_body = (entity.radius * 2.0 * self.grid_ratio).ceil() as usize + 1;
        entity.grid_pos_x = ((entity.x - entity.radius) * self.grid_ratio) as usize;
        entity.grid_pos_y = ((entity.y - entity.radius) * self.grid_ratio) as usize;

        for y_pos in entity.grid_pos_y
            ..usize::min(entity.grid_pos_y + entity.grid_body, self.spatial_grid_dimension)
        {
            for x_pos in entity.grid_pos_x
                ..usize::min(entity.grid_pos_x + entity.grid_body, self.spatial_grid_dimension)
            {
                update_grid_position(
                    (y_pos << self.encoding_bits) | x_pos,
                    index,
                    spatial_grid_ptr,
                    spatial_grid_locks_ptr,
                );
            }
        }
        entity.radius
    }

    #[inline(always)]
    pub fn remove_entity(&mut self, index: usize) {
        let entity = &mut self.entities[index];
        if entity.replace {
            return;
        };
        entity.replace = true;

        self.replacement_queue.push(entity.index);

        let spatial_grid_ptr = self.spatial_grid.as_mut_ptr();
        let spatial_grid_locks_ptr = self.spatial_grid_locks.as_mut_ptr();

        for y_pos in entity.grid_pos_y
            ..usize::min(entity.grid_pos_y + entity.grid_body, self.spatial_grid_dimension)
        {
            for x_pos in entity.grid_pos_x
                ..usize::min(entity.grid_pos_x + entity.grid_body, self.spatial_grid_dimension)
            {
                remove_grid_position(
                    (y_pos << self.encoding_bits) | x_pos,
                    &entity.index,
                    spatial_grid_ptr,
                    spatial_grid_locks_ptr,
                );
            }
        }
    }

    #[inline(always)]
    pub fn create_entity_movement_from_angle(&mut self, index: usize, angle: f32) {
        let entity = &mut self.entities[index];
        entity.acceleration_x = f32::cos(angle) * entity.movement_acceleration;
        entity.acceleration_y = f32::sin(angle) * entity.movement_acceleration;
    }

    #[inline(always)]
    pub fn create_entity_movement_from_cardinal_direction(&mut self, index: usize, direction: u8) {
        let entity = &mut self.entities[index];
        match direction {
            0 => {
                entity.acceleration_x = entity.movement_acceleration;
                entity.acceleration_y = 0.0;
            }
            1 => {
                entity.acceleration_x = -entity.movement_acceleration;
                entity.acceleration_y = 0.0;
            }
            2 => {
                entity.acceleration_y = entity.movement_acceleration;
                entity.acceleration_x = 0.0;
            }
            3 => {
                entity.acceleration_y = -entity.movement_acceleration;
                entity.acceleration_x = 0.0;
            }
            4 => {
                entity.acceleration_x = entity.movement_acceleration * 0.70710675;
                entity.acceleration_y = entity.movement_acceleration * 0.70710675;
            }
            5 => {
                entity.acceleration_x = entity.movement_acceleration * 0.70710675;
                entity.acceleration_y = -entity.movement_acceleration * 0.70710675;
            }
            6 => {
                entity.acceleration_x = -entity.movement_acceleration * 0.70710675;
                entity.acceleration_y = entity.movement_acceleration * 0.70710675;
            }
            7 => {
                entity.acceleration_x = -entity.movement_acceleration * 0.70710675;
                entity.acceleration_y = -entity.movement_acceleration * 0.70710675;
            }
            _ => {
                entity.acceleration_x = 0.0;
                entity.acceleration_y = 0.0;
            }
        }
    }

    #[inline(always)]
    pub fn update(&mut self) {
        let room_size = self.size;
        let room_grid_ratio = self.grid_ratio;
        let entities_size = self.entities.len();

        let entities_ptr = self.entities.as_mut_ptr() as usize;
        let spatial_grid_ptr = self.spatial_grid.as_mut_ptr() as usize;
        let spatial_grid_locks_ptr = self.spatial_grid_locks.as_mut_ptr() as usize;

        let spatial_grid_dimension = self.spatial_grid_dimension;
        let encoding_bits = self.encoding_bits;

        rayon::scope(|s| {
            for thread in 0..THREADS {
                let chunk = self.entity_chunks[thread];
                s.spawn(move |_| {
                    let entities = entities_ptr as *mut Entity;
                    let spatial_grid = spatial_grid_ptr as *mut HashSet<usize>;
                    let spatial_grid_locks = spatial_grid_locks_ptr as *mut AtomicBool;
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
                            } else if entity.x > room_size - entity.radius {
                                entity.x = room_size - entity.radius;
                                entity.velocity_x = -entity.velocity_x;
                            };
                            if entity.y < entity.radius {
                                entity.y = entity.radius;
                                entity.velocity_y = -entity.velocity_y;
                            } else if entity.y > room_size - entity.radius {
                                entity.y = room_size - entity.radius;
                                entity.velocity_y = -entity.velocity_y;
                            };
                            entity.velocity_x *= entity.friction;
                            entity.velocity_y *= entity.friction;
                            let speed_sq = entity.velocity_x * entity.velocity_x
                                + entity.velocity_y * entity.velocity_y;
                            if speed_sq
                                > entity.terminal_velocity_in_direction
                                    * entity.terminal_velocity_in_direction
                            {
                                let speed_scalar =
                                    entity.terminal_velocity_in_direction / speed_sq.sqrt();
                                entity.velocity_x *= speed_scalar;
                                entity.velocity_y *= speed_scalar;
                            };
                            let spatial_grid_pos_x =
                                ((entity.x - entity.radius) * room_grid_ratio) as usize;
                            let spatial_grid_pos_y =
                                ((entity.y - entity.radius) * room_grid_ratio) as usize;
                            if spatial_grid_pos_x != entity.grid_pos_x
                                || spatial_grid_pos_y != entity.grid_pos_y
                            {
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
                                    update_entity_body(
                                        entity.grid_body,
                                        entity.grid_body,
                                        spatial_grid_pos_y,
                                        spatial_grid_pos_x,
                                        entity.grid_pos_y,
                                        entity.grid_pos_x,
                                        entity.index,
                                        spatial_grid,
                                        spatial_grid_locks,
                                        spatial_grid_dimension,
                                        encoding_bits,
                                    );
                                } else {
                                    if shift_x.0 == 0 {
                                        if shift_y.1 {
                                            update_entity_body(
                                                shift_y.0,
                                                entity.grid_body,
                                                entity.grid_pos_y + entity.grid_body,
                                                spatial_grid_pos_x,
                                                entity.grid_pos_y,
                                                spatial_grid_pos_x,
                                                entity.index,
                                                spatial_grid,
                                                spatial_grid_locks,
                                                spatial_grid_dimension,
                                                encoding_bits,
                                            );
                                        } else {
                                            update_entity_body(
                                                shift_y.0,
                                                entity.grid_body,
                                                spatial_grid_pos_y,
                                                spatial_grid_pos_x,
                                                spatial_grid_pos_y + entity.grid_body,
                                                spatial_grid_pos_x,
                                                entity.index,
                                                spatial_grid,
                                                spatial_grid_locks,
                                                spatial_grid_dimension,
                                                encoding_bits,
                                            );
                                        };
                                    } else if shift_y.0 == 0 {
                                        if shift_x.1 {
                                            update_entity_body(
                                                entity.grid_body,
                                                shift_x.0,
                                                spatial_grid_pos_y,
                                                entity.grid_pos_x + entity.grid_body,
                                                spatial_grid_pos_y,
                                                entity.grid_pos_x,
                                                entity.index,
                                                spatial_grid,
                                                spatial_grid_locks,
                                                spatial_grid_dimension,
                                                encoding_bits,
                                            );
                                        } else {
                                            update_entity_body(
                                                entity.grid_body,
                                                shift_x.0,
                                                spatial_grid_pos_y,
                                                spatial_grid_pos_x,
                                                spatial_grid_pos_y,
                                                spatial_grid_pos_x + entity.grid_body,
                                                entity.index,
                                                spatial_grid,
                                                spatial_grid_locks,
                                                spatial_grid_dimension,
                                                encoding_bits,
                                            );
                                        };
                                    } else {
                                        if shift_y.1 {
                                            update_entity_body(
                                                shift_y.0,
                                                entity.grid_body,
                                                spatial_grid_pos_y + entity.grid_body - shift_y.0,
                                                spatial_grid_pos_x,
                                                entity.grid_pos_y,
                                                entity.grid_pos_x,
                                                entity.index,
                                                spatial_grid,
                                                spatial_grid_locks,
                                                spatial_grid_dimension,
                                                encoding_bits,
                                            );
                                        } else {
                                            update_entity_body(
                                                shift_y.0,
                                                entity.grid_body,
                                                spatial_grid_pos_y,
                                                spatial_grid_pos_x,
                                                entity.grid_pos_y + entity.grid_body - shift_y.0,
                                                entity.grid_pos_x,
                                                entity.index,
                                                spatial_grid,
                                                spatial_grid_locks,
                                                spatial_grid_dimension,
                                                encoding_bits,
                                            );
                                        };
                                        let y_overlap_start =
                                            usize::max(spatial_grid_pos_y, entity.grid_pos_y);
                                        if shift_x.1 {
                                            update_entity_body(
                                                entity.grid_body - shift_y.0,
                                                shift_x.0,
                                                y_overlap_start,
                                                spatial_grid_pos_x + entity.grid_body - shift_x.0,
                                                y_overlap_start,
                                                entity.grid_pos_x,
                                                entity.index,
                                                spatial_grid,
                                                spatial_grid_locks,
                                                spatial_grid_dimension,
                                                encoding_bits,
                                            );
                                        } else {
                                            update_entity_body(
                                                entity.grid_body - shift_y.0,
                                                shift_x.0,
                                                y_overlap_start,
                                                spatial_grid_pos_x,
                                                y_overlap_start,
                                                entity.grid_pos_x + entity.grid_body - shift_x.0,
                                                entity.index,
                                                spatial_grid,
                                                spatial_grid_locks,
                                                spatial_grid_dimension,
                                                encoding_bits,
                                            );
                                        };
                                    };
                                };
                                entity.grid_pos_x = spatial_grid_pos_x;
                                entity.grid_pos_y = spatial_grid_pos_y;
                            };
                        };
                    }
                });
            }
        });

        let required_len = THREADS * entities_size;
        if self.applied_collisions.capacity() < required_len {
            self.applied_collisions
                .reserve(required_len - self.applied_collisions.len());
        }
        unsafe { self.applied_collisions.set_len(required_len) };
        let applied_collisions_ptr = self.applied_collisions.as_mut_ptr() as usize;

        rayon::scope(|s| {
            for thread in 0..THREADS {
                s.spawn(move |_| {
                    let applied_collisions = applied_collisions_ptr as *mut [f32; 4];
                    unsafe {
                        let bucket_start = applied_collisions.add(thread * entities_size);
                        ptr::write_bytes(bucket_start, 0, entities_size);
                    };
                });
            }
        });

        let stored_collisions_ptr = self.stored_collisions.as_mut_ptr() as usize;

        rayon::scope(|s| {
            for thread in 0..THREADS {
                let chunk = self.entity_chunks[thread];
                s.spawn(move |_| {
                    let applied_collisions = applied_collisions_ptr as *mut [f32; 4];
                    let entities = entities_ptr as *mut Entity;
                    let spatial_grid = spatial_grid_ptr as *mut HashSet<usize>;

                    let stored_collisions_slice = stored_collisions_ptr as *mut Vec<[usize; 2]>;
                    let stored_collisions = unsafe { &mut *stored_collisions_slice.add(thread) };
                    stored_collisions.clear();

                    unsafe {
                        for entity_index in chunk[0]..chunk[1] {
                            let entity = &*entities.add(entity_index);
                            if entity.replace {
                                continue;
                            };
                            for y_pos in entity.grid_pos_y
                                ..usize::min(
                                    entity.grid_pos_y + entity.grid_body,
                                    spatial_grid_dimension,
                                )
                            {
                                for x_pos in entity.grid_pos_x
                                    ..usize::min(
                                        entity.grid_pos_x + entity.grid_body,
                                        spatial_grid_dimension,
                                    )
                                {
                                    let grid_index = (y_pos << encoding_bits) | x_pos;
                                    let spatial_grid_position = &*spatial_grid.add(grid_index);
                                    if spatial_grid_position.len() <= 1 {
                                        continue;
                                    };

                                    for &collision_index in spatial_grid_position.iter() {
                                        if entity_index >= collision_index {
                                            continue;
                                        };
                                        let collision_entity = &*entities.add(collision_index);
                                        let min_shared_x = usize::max(
                                            entity.grid_pos_x,
                                            collision_entity.grid_pos_x,
                                        );
                                        let min_shared_y = usize::max(
                                            entity.grid_pos_y,
                                            collision_entity.grid_pos_y,
                                        );
                                        let min_shared_coordinate =
                                            (min_shared_y << encoding_bits) | min_shared_x;

                                        if min_shared_coordinate == grid_index {
                                            if entity.body_type + collision_entity.body_type == 2 {
                                                let dsq = (collision_entity.x - entity.x)
                                                    * (collision_entity.x - entity.x)
                                                    + (collision_entity.y - entity.y)
                                                        * (collision_entity.y - entity.y);
                                                if dsq
                                                    < (entity.radius + collision_entity.radius)
                                                        * (entity.radius + collision_entity.radius)
                                                {
                                                    manage_collision_circle_circle(
                                                        dsq,
                                                        entity,
                                                        collision_entity,
                                                        applied_collisions,
                                                        thread,
                                                        entities_size,
                                                        stored_collisions,
                                                    );
                                                };
                                            } else if entity.body_type + collision_entity.body_type
                                                == 1
                                            {
                                                if entity.body_type == 0 {
                                                    let cx = collision_entity.x.clamp(
                                                        entity.x - entity.radius,
                                                        entity.x + entity.radius,
                                                    );
                                                    let cy = collision_entity.y.clamp(
                                                        entity.y - entity.radius,
                                                        entity.y + entity.radius,
                                                    );
                                                    if (cx - collision_entity.x)
                                                        * (cx - collision_entity.x)
                                                        + (cy - collision_entity.y)
                                                            * (cy - collision_entity.y)
                                                        < collision_entity.radius
                                                            * collision_entity.radius
                                                    {
                                                        manage_collision_circle_square(
                                                            cx,
                                                            cy,
                                                            collision_entity,
                                                            entity,
                                                            applied_collisions,
                                                            thread,
                                                            entities_size,
                                                            stored_collisions,
                                                        );
                                                    };
                                                } else {
                                                    let cx = entity.x.clamp(
                                                        collision_entity.x
                                                            - collision_entity.radius,
                                                        collision_entity.x
                                                            + collision_entity.radius,
                                                    );
                                                    let cy = entity.y.clamp(
                                                        collision_entity.y
                                                            - collision_entity.radius,
                                                        collision_entity.y
                                                            + collision_entity.radius,
                                                    );
                                                    if (cx - entity.x) * (cx - entity.x)
                                                        + (cy - entity.y) * (cy - entity.y)
                                                        < entity.radius * entity.radius
                                                    {
                                                        manage_collision_circle_square(
                                                            cx,
                                                            cy,
                                                            entity,
                                                            collision_entity,
                                                            applied_collisions,
                                                            thread,
                                                            entities_size,
                                                            stored_collisions,
                                                        );
                                                    };
                                                };
                                            } else {
                                                if f32::abs(collision_entity.x - entity.x)
                                                    < entity.radius + collision_entity.radius
                                                    && f32::abs(collision_entity.y - entity.y)
                                                        < entity.radius + collision_entity.radius
                                                {
                                                    manage_collision_square_square(
                                                        entity,
                                                        collision_entity,
                                                        applied_collisions,
                                                        thread,
                                                        entities_size,
                                                        stored_collisions,
                                                    );
                                                };
                                            };
                                        }
                                    }
                                }
                            }
                        }
                    }
                });
            }
        });

        rayon::scope(|s| {
            for thread in 0..THREADS {
                let chunk = self.entity_chunks[thread];
                s.spawn(move |_| {
                    let applied_collisions = applied_collisions_ptr as *mut [f32; 4];
                    let entities = entities_ptr as *mut Entity;
                    unsafe {
                        for index in chunk[0]..chunk[1] {
                            let entity = &mut *entities.add(index);
                            if entity.replace {
                                continue;
                            };
                            let mut total_simd = f32x4::ZERO;
                            for bucket in 0..THREADS {
                                let bucket_pos =
                                    applied_collisions.add(bucket * entities_size + index);
                                total_simd += f32x4::from(*bucket_pos);
                            }
                            let total = total_simd.to_array();
                            entity.velocity_x += total[0];
                            entity.velocity_y += total[1];
                            entity.x += total[2];
                            entity.y += total[3];
                        }
                    };
                });
            }
        });

        self.tick += 1;
    }
}
