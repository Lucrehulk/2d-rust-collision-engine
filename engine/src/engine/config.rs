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
