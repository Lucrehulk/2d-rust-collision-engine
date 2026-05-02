use collision_engine::engine::engine::Room;
use std::time::{Instant, Duration};
use rand::Rng;

// --- Benchmark Configuration Constants ---
const NUM_CIRCLES: usize = 5000;
const NUM_SQUARES: usize = 5000;
const ROOM_SIZE: f32 = 1024.0;

fn main() {
    let mut world = Room::init(ROOM_SIZE);
    let mut rng = rand::thread_rng();

    println!("Spawning {} circles and {} squares...", NUM_CIRCLES, NUM_SQUARES);

    // Spawn Circles (body_type = 0)
    for _ in 0..NUM_CIRCLES {
        let radius = rng.gen_range(2.0..=6.0);
        let x = rng.gen_range(radius..(ROOM_SIZE - radius));
        let y = rng.gen_range(radius..(ROOM_SIZE - radius));
        let vx = rng.gen_range(-2.0..=2.0);
        let vy = rng.gen_range(-2.0..=2.0);
        let terminal_velocity = 3.0; 
        
        let id = world.create_entity(
            x, y, 
            10.0, // mass
            0.99, // friction
            vx, vy, 
            terminal_velocity, 
            0.1,  // movement_acceleration
            radius, 
            0,    // body_type: 0 for circle
            1     // collision_type
        ); //

        // Manually assign a random acceleration since create_entity defaults it to 0.0[cite: 3]
        world.entities[id].acceleration_x = rng.gen_range(-0.05..=0.05);
        world.entities[id].acceleration_y = rng.gen_range(-0.05..=0.05);
    }

    // Spawn Squares (body_type = 1)
    for _ in 0..NUM_SQUARES {
        let radius = rng.gen_range(2.0..=6.0);
        let x = rng.gen_range(radius..(ROOM_SIZE - radius));
        let y = rng.gen_range(radius..(ROOM_SIZE - radius));
        let vx = rng.gen_range(-2.0..=2.0);
        let vy = rng.gen_range(-2.0..=2.0);
        let terminal_velocity = 3.0; 
        
        let id = world.create_entity(
            x, y, 
            10.0, // mass
            0.99, // friction
            vx, vy, 
            terminal_velocity, 
            0.1,  // movement_acceleration
            radius, 
            1,    // body_type: 1 for square
            1     // collision_type
        ); //[cite: 3]

        world.entities[id].acceleration_x = rng.gen_range(-0.05..=0.05);
        world.entities[id].acceleration_y = rng.gen_range(-0.05..=0.05);
    }

    println!("Initialization complete. Starting benchmark...");

    let mut total_duration = Duration::new(0, 0);
    let mut tick_count = 0;

    // Run the benchmark continuously
    loop {
        let start = Instant::now();
        
        world.update(); //[cite: 3]
        
        let elapsed = start.elapsed();
        total_duration += elapsed;
        tick_count += 1;
        
        let average_duration = total_duration / tick_count;

        println!(
            "Tick: {:0>5} | Last Tick Time: {:>8.2?} | Average Time: {:>8.2?}",
            tick_count, elapsed, average_duration
        );
    }
}