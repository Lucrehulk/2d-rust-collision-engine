# 2d-rust-collision-engine
A simple and fast 2d collision managing, physics system system, designed to be utilized for io games.

Note that this system does not support collision management for all polygon types--only circles and axis-aligned squares (for simple hitbox detection).

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
    // Determines if friction is enabled.
    pub const FRICTION: bool = false;
    // The constant of friction, this is consistently multiplied to the entity's x and y velocity every tick.
    pub const FRICTIONAL_CONSTANT: f32 = 0.99;
    // Determines if gravity is enabled.
    pub const GRAVITY: bool = false;
    // The constant of the acceleration due to gravity, this is consistently added to the entity's y velocity so as to make the entity go downwards (hence replicating gravity).
    pub const GRAVITY_CONSTANT: f32 = 0.2;
    // The acceleration applied to both the x and y velocities when a collision occurs. 
    pub const COLLISION_ACCELERATION_CONSTANT: f32 = 0.2;
    // The acceleration applied to an entity when it is in a movement state.
    pub const MOVEMENT_ACCELERATION_CONSTANT: f32 = 0.01;
    // The downtime, in ms, per tick.
    pub const TICK_TIME: u64 = 16;
}
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
50000 entities: ~15mspt
100000 entities: ~30-40mspt (with so many entities this can often vary between ticks by a few)
```

Additionally, other features such as spawning/deleting entities, support for creating self-entity driven movement (via cardinal directions or given a specific angle), and more are included.
