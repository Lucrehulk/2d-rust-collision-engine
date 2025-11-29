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
