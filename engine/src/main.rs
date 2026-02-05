// Demo server. 
use collision_engine::engine::engine::Room;
use collision_engine::engine::config::config_data::{
    ROOM_SIZE,
    SPATIAL_GRID_DIMENSION,
    TICK_TIME,
    THREADS,
    STORE_COLLISIONS,
};
use rand::Rng; 
use std::{time::{Duration, Instant}, sync::Arc};
use tokio::{net::TcpListener, sync::broadcast, time::sleep};
use tokio_tungstenite::{accept_async, tungstenite::Message};
use futures_util::{StreamExt, SinkExt};

#[tokio::main]
async fn main() {
    let mut world = Room::init();
    let mut rng = rand::rng();

    let shape_group_count = 5000;
    
    for _ in 0..shape_group_count {
        // Generate random square body
        let x = rng.random_range(0.0..ROOM_SIZE);
        let y = rng.random_range(0.0..ROOM_SIZE);
        let mass = rng.random_range(1.0..2.0);
        let vx = rng.random_range(-2.0..2.0); 
        let vy = rng.random_range(-2.0..2.0); 
        let radius = rng.random_range(2.0..6.0);
        world.create_entity(x, y, mass, 0.0, vx, vy, 2.0, 0.0, radius, 0, false);
        // Generate random circle body
        let x = rng.random_range(0.0..ROOM_SIZE);
        let y = rng.random_range(0.0..ROOM_SIZE);
        let mass = rng.random_range(1.0..2.0);
        let vx = rng.random_range(-2.0..2.0); 
        let vy = rng.random_range(-2.0..2.0); 
        let radius = rng.random_range(2.0..6.0);
        world.create_entity(x, y, mass, 0.0, vx, vy, 2.0, 0.0, radius, 1, false);
    }

    // world.remove_entity(0);
    // world.remove_entity(1);
    // world.create_entity(ROOM_SIZE / 2.0, ROOM_SIZE / 2.0, f32::MAX, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0, false);

    let listener = TcpListener::bind("127.0.0.1:8080").await.unwrap();
    println!("Server running at ws://127.0.0.1:8080");
    let (tx, _) = broadcast::channel::<Vec<u8>>(16);
    let tx = Arc::new(tx);
    let tx_clone = tx.clone();
    tokio::spawn(async move {
        loop {
            match listener.accept().await {
                Ok((stream, _)) => {
                    let tx = tx_clone.clone();
                    tokio::spawn(async move {
                        let ws_stream = match accept_async(stream).await {
                            Ok(ws) => ws,
                            Err(e) => {
                                eprintln!("WebSocket error: {}", e);
                                return;
                            }
                        };

                        let (mut write, _) = ws_stream.split();
                        let mut rx = tx.subscribe();

                        let mut packet = vec![0]; 
                        packet.extend(ROOM_SIZE.to_le_bytes());
                        packet.extend((SPATIAL_GRID_DIMENSION as u32).to_le_bytes());
                        if write.send(Message::Binary(packet)).await.is_err() {
                            return;
                        }

                        while let Ok(packet) = rx.recv().await {
                            if write.send(Message::Binary(packet.clone())).await.is_err() {
                                break;
                            }
                        }
                    });
                }
                Err(e) => eprintln!("Failed to accept connection: {}", e),
            }
        }
    });

    let mut tick_sum = 0.0;

    loop {
        let tick = Instant::now();
        world.update();
        let tick_time = tick.elapsed();
        tick_sum += tick_time.as_secs_f64() * 1000.0;
        println!("tick #: {:?}, tick time: {:?}, avg tick time: {:?}", world.tick, tick_time, tick_sum / (world.tick as f64));

        // Set STORE_COLLISIONS in config to true to see stored collision data. This can be used for external game logic if necessary.
        if STORE_COLLISIONS {
            for thread in 0..THREADS {
                println!("{} - {}", thread, world.stored_collisions[thread].len());
                for collision in &world.stored_collisions[thread] {
                    println!("this is just a test showing there was a collision between entities {} and {} when STORE_COLLISIONS config value was set to true.", collision[0], collision[1]);
                }
            }
        }

        // Packet stuff, pretty unoptimized as it just sends all entity data but this is just a proof of concept.
        let mut packet = vec![1];
        for entity in &world.entities {
            if !entity.replace {
                packet.extend((entity.index as u32).to_le_bytes());
                packet.extend(entity.x.to_le_bytes());
                packet.extend(entity.y.to_le_bytes());
                packet.extend(entity.radius.to_le_bytes());
                packet.push(entity.body_type);
            }
        }
        let _ = tx.send(packet);
        
        sleep(Duration::from_millis(TICK_TIME)).await;
    }
}
