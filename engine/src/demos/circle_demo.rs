mod config;
mod engine;
use crate::engine::Room;
use crate::config::config_data::{
    ROOM_SIZE,
    SPATIAL_GRID_DIMENSION,
    THREADS,
    MAX_ENTITIES,
    MAX_ENTITIES_TO_REPLACE,
    FRICTION,
    FRICTIONAL_CONSTANT,
    GRAVITY,
    GRAVITY_CONSTANT,
    COLLISION_ACCELERATION_CONSTANT,
    MOVEMENT_ACCELERATION_CONSTANT,
    TICK_TIME
};
use rand::Rng; 
use std::{thread, collections::HashSet, time::{Duration, Instant}, sync::Mutex};

// Just example usage and stuff
fn main() {
    let mut world = Room::init();
    world.collision_positions_ptr = &mut world.collision_positions as *mut Mutex<HashSet<usize>>;
    
    // Create random circle entities example
    let mut rng = rand::rng();

    for _ in 0..MAX_ENTITIES {
        let x = rng.random_range(0.0..ROOM_SIZE);
        let y = rng.random_range(0.0..ROOM_SIZE);
        let vx = rng.random_range(-2.0..2.0); 
        let vy = rng.random_range(-2.0..2.0); 
        let radius = rng.random_range(2.0..6.0);
        let body_type = 1;
        world.create_entity(x, y, vx, vy, 2.0, 2.0, radius, body_type);
    }
    
    loop {
        let tick = Instant::now();
        world.update();
        let tick_time = tick.elapsed();
        println!("{:?}", tick_time);
        thread::sleep(Duration::from_millis(TICK_TIME));
    }
}
