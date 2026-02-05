# 2d-rust-collision-engine
A simple and fast 2d collision managing, physics elastic collision system, designed to be utilized for io games.

Note that this system does not support collision management for all polygon types--only circles and axis-aligned squares (for simple hitbox detection).

The engine is primarily designed to handle large amounts of entities, and uses multithreading techniques in order to effectively do so. For environments with few entities, it is rather slow due to this. The engine performs best in medium sized enviroments that contain large amounts of entities.

For those wanting specific numbers, I did benchmark the performance of the system in an environment with these config settings:
```
pub mod config_data {
    // Room size.
    pub const ROOM_SIZE: f32 = 1024.0;
    // Amount of spatial grid units on a single side (i.e. area = SPATIAL_GRID_DIMENSION * SPATIAL_GRID_DIMENSION).
    pub const SPATIAL_GRID_DIMENSION: usize = 64;
    // Target number of threads for the application to use.
    pub const THREADS: usize = 16;
    // An initial capacity for the max number of entities the system is designed to handle. This value can be exceeded but if you know the maximum amount of entities in your system and are willing to allocate the space it is slightly more optimal to set a known capacity.
    pub const MAX_ENTITIES: usize = 2;
    // The capacity for the maximum number of open entity indexes the system is designed to handle. This is just like the MAX_ENTITIES constant, it can be exceeded and is not necessary. It is a little optimization if you know the general maximum value you'll be working with though, and if you'll be around that value often.
    pub const MAX_ENTITIES_TO_REPLACE: usize = 2;
    // The downtime, in ms, per tick.
    pub const TICK_TIME: u64 = 16;
    // The slight separation constant added when dealing with displacements for collisions. This is added so there is no risk entities still will just be barely touching. 
    pub const COLLISION_SEPARATION_CONSTANT: f32 = 0.00001;
    // Option to store the ids of entities in collisions within the stored_collisions field of Room. Useful if you need collisions for logic outside of physics (eg. hp in a game).
    pub const STORE_COLLISIONS: bool = false;
}
```

The device this system was benchmarked on is as follows:
```
Processor	11th Gen Intel(R) Core(TM) i7-11700 @ 2.50GHz   2.50 GHz
Installed RAM	16.0 GB (15.7 GB usable)
System type	64-bit operating system, x64-based processor
Pen and touch	No pen or touch input is available for this display
```

Field data was generated as so for all tests:
```
let x = rng.random_range(0.0..ROOM_SIZE);
let y = rng.random_range(0.0..ROOM_SIZE);
let mass = rng.random_range(1.0..2.0);
let vx = rng.random_range(-2.0..2.0); 
let vy = rng.random_range(-2.0..2.0); 
let radius = rng.random_range(2.0..6.0);

terminal_velocity_in_direction was always 2.0.
Friction disabled for all entities.
```

ENTITIES 1/2 CIRCLE, 1/2 SQUARE MIXED BENCHMARK
```
100 entities -> 2.9216670000000007
1000 entities -> 3.322058
10000 entities -> 4.709827000000001
20000 entities -> 7.795267999999998
50000 entities -> 25.044604000000003
100000 entities -> 86.77300499999996
```

ENTITIES ALL CIRCLE BENCHMARK
```
100 entities -> 3.245687999999999
1000 entities -> 3.5085800000000003
10000 entities -> 4.945434999999998
20000 entities -> 7.323725000000003
50000 entities -> 21.677546999999993
100000 entities -> 70.83480400000002
```

ENTITIES ALL SQUARE BENCHMARK
```
100 entities -> 3.0125319999999998
1000 entities -> 3.552559000000001
10000 entities -> 5.004377000000001
20000 entities -> 7.649921000000001
50000 entities -> 23.304519
100000 entities -> 79.928201
```

For your own testing, a basic demo server (main.rs) and html client are provided. You can view and benchmark the system in action with these.
