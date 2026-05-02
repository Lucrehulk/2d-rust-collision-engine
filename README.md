# 2d-rust-collision-engine
A simple and fast 2d collision managing, physics elastic collision system, designed to be utilized for io games.

This is my own personal physics engine I use.

Note that this system does not support collision management for all polygon types--only circles and axis-aligned squares (for simple hitbox detection).

The engine is primarily designed to handle large amounts of entities, and uses multithreading techniques in order to effectively do so.

Engine utilizes:

```
Spatial hash grid implementation.

Impulse collision calculations for minimal overhead.

Threads and chunk division of work to avoid most slower locking structures. Instead, data separation/isolation, along with atomics to avoid traditional locking. Instead, we simply use our own custom "locking" via atomic booleans, which determine if the grid position is occupied or not.

Thread pooling to minimize the overhead of constantly starting and stopping threads.

Resizable vector for entities, with ID mapping to index and replaceable ID queue to avoid HashMap usage.

Partial grid body updates to avoid having to update unnecessary entity grid body positions in the spatial grid.

Various unsafe methods for minor speedup and optimization.
```



I have not yet thoroughly benchmarked the engine's performance, but it works well. A benchmark test file is provided. I ran 5k circles and 5k squares on my i7 8 core device with 16 threads and had it running under 6ms per tick. 

Overall, can it be used for games? Yes. It's also pretty fast and much simpler compared to most engines. It's not the fastest engine ever, but it can still handle tens of thousands of entities in reasonable time.

Further updates like more advanced collision ignoring typing (ability to ignore or enact physical collisions based on a given type ID), and potential optimizations like thread pooling may also be added soon.
