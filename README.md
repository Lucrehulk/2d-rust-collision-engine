# 2d-rust-collision-engine
A simple and fast 2d collision managing, physics system system, designed to be utilized for io games.

Note that this system does not support collision management for all polygon types--only circles and axis-aligned squares (for simple hitbox detection).

At the top of the collision system's code, constants for configuration (i.e. room data, physics config like gravity, friction, and collision accelerations) can be found.

The demo server and client already provides a basic demo of the system in action, while the physics_engine.rs file provides the raw, non-networked code for the physics engine only. 

Some examples from the demo system:
10k entities in a 1024x1024 space, with a giant immovable wall in the middle:
![collision demo](https://github.com/user-attachments/assets/56394174-3f88-455f-ba57-0d3d2f20125b)
A video with 1k entities (for easier viewing of the physics) in the same 1024x1024 space, with the same immovable wall:
[collision demo video](https://github.com/user-attachments/assets/c5aac9a9-3a0b-4179-9632-15f5436b3144)

For those wanting specific numbers, I did benchmark the performance of the system in an environment with these settings:
```
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
```

The device this system was benchmarked on is as follows:
```
Processor	11th Gen Intel(R) Core(TM) i7-11700 @ 2.50GHz   2.50 GHz
Installed RAM	16.0 GB (15.7 GB usable)
System type	64-bit operating system, x64-based processor
Pen and touch	No pen or touch input is available for this display
```

The entities themselves were circular only (though square entities would yield similar performance still), with their x and y positions being randomly generated to be anywhere within the room, x and y velocities were within -2 and 2 (maximum magnitude of velocity for every entity being 2), and the radius of each entity was in between 2 and 6.

The results were as follows:
```
100 entities: ~1.5mspt
1000 entities: ~2mspt
10000 entities: ~3mspt
Note that at this point, the 1024x1024 space can't even effectively fit all the entities without some overlap anymore, an image of the 50000 entities is provided below. Note that this is why such a drastic increase in mspt occurs here. In a larger space (as I will show in another benchmark, the performance is still silky smooth)
50000 entities: ~10-20mspt
100000 entities: ~45-50mspt
```

An example of how cluttered the system becomes at only 50k entities in the 1024x1024 environment:
![Screenshot 2025-06-28 144716](https://github.com/user-attachments/assets/85299206-d514-4d65-8549-34ac7275c7d8)

However, for example, after extending the room size to 4096x4096, and extending the dimensions of the spatial grid to 128x128, I found that 100000 entities could perform at even under 10mspt for certain ticks. Note that the environment and spatial grid is still rather cluttered with that many entities.

A snippet of the 100k entities in the 4096x4096 room, 128x128 spatial grid benchmark:
![Screenshot 2025-06-28 145522](https://github.com/user-attachments/assets/97f36325-de1a-4120-88cb-5aac630147eb)
