# 2d-rust-collision-engine
A simple and fast 2d collision managing, physics elastic collision system, designed to be utilized for io games.

This is my own personal physics engine I use.

Note that this system does not support collision management for all polygon types--only circles and axis-aligned squares (for simple hitbox detection).

The engine is primarily designed to handle large amounts of entities, and uses multithreading techniques in order to effectively do so.

Engine utilizes:

---

Standard position update phase -> collision detection phase -> force/impulse updates phase (update position & velocity from collisions).

Spatial hash grid implementation (2d vector of hashsets containing ids).

Impulse collision calculations for minimal overhead.

Threads and chunk division of work to avoid most slower locking structures. Instead, data separation/isolation. Entities vector is chunked into as even chunks as possible for the threads to perform work on individually. This is used across the three core parts of the engine. 

Atomics to avoid traditional locking. Instead, we simply use our own custom "locking" via atomic booleans, which determine if the grid position is occupied or not (a sizeXsize spatial grid map dimension copy for a 2d array of atomic booleans is implemented. When adding or removing an id from the spatial grid, we do a `compare_exchange_weak(false, true` and `spin_loop()` while this fails. exchanging the false to true lets the rest of the threads reading at that spatial grid know that index is currently being accessed, and thise other threads will keep spin looping until the thread currently holding the spatial grid location uses `store` to set the atomic boolean back to false. It is effectively a much faster locking method.

Thread pooling to minimize the overhead of constantly starting and stopping threads. Rayon.

Resizable vector for entities, with ID mapping to index and replaceable ID queue to avoid HashMap usage. Effectively entities are a nice vector, with them taking ids based on vector index. Deleted entities are not actually deleted but instead given a flag marking them for deletion, so they don't get updated or checked for collisions. The vector for the queue, well, queues these indices and eventually replaces them with new entities assuming more entity spawn calls are made.

Partial grid body updates to avoid having to update unnecessary entity grid body positions in the spatial grid. Effectively, the bodies are treated as squares (since they are radius/1 side based, no rectangles) in the spatial grid. Using some basic math involving the ahifts in these grid bodies, its checked if an update is needed to the spatial grid at all (coordinates changed, entity resize method automatically updates it). If so, iteration bounds that calculate indices to iterate the two non-overlapping regions (deletion region + update region) are calculated by using the original position and shift math along with the body size. 

Various unsafe methods for minor speedup and optimization. Would've just used C++ but that's not on my stack unfortunately.

---

I have not yet thoroughly benchmarked the engine's performance, but it works well. A benchmark test file is provided. I ran 5k circles and 5k squares on my i7 8 core device with 16 threads and had it running under 6ms per tick. 

Overall, can it be used for games? Yes. It's also pretty fast and much simpler compared to most engines. It's not the fastest engine ever, but it can still handle tens of thousands of entities in reasonable time.
