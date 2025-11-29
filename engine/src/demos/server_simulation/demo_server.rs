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
use std::{thread, collections::HashSet, time::{Duration, Instant}, sync::{Mutex, Arc}};
use tokio::{net::TcpListener, sync::broadcast, time::sleep};
use tokio_tungstenite::{accept_async, tungstenite::Message};
use futures_util::{StreamExt, SinkExt};

#[tokio::main]
async fn main() {
    let mut world = Room::init();
    world.collision_positions_ptr = &mut world.collision_positions as *mut Mutex<HashSet<usize>>;
    let mut rng = rand::rng();

    for _ in 0..MAX_ENTITIES {
        let x = rng.random_range(0.0..ROOM_SIZE);
        let y = rng.random_range(0.0..ROOM_SIZE);
        let vx = rng.random_range(-2.0..2.0); 
        let vy = rng.random_range(-2.0..2.0); 
        let radius = rng.random_range(2.0..6.0);
        let body_type = f32::round(rng.random_range(0.0..1.0)) as u8;
        world.create_entity(x, y, vx, vy, 2.0, 2.0, radius, body_type);
    }
    world.remove_entity(0);
    world.create_entity(ROOM_SIZE / 2.0, ROOM_SIZE / 2.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0);

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

    loop {
        let tick = Instant::now();
        world.update();
        let tick_time = tick.elapsed();
        println!("{:?}", tick_time);

        // Packet stuff, pretty unoptimized as it just sends all entity data but this is just a proof of concept.
        let mut packet = vec![1];
        for entity in &world.entities {
            packet.extend((entity.index as u32).to_le_bytes());
            packet.extend(entity.x.to_le_bytes());
            packet.extend(entity.y.to_le_bytes());
            packet.extend(entity.radius.to_le_bytes());
            packet.push(entity.body_type);
        }
        let _ = tx.send(packet);
        
        sleep(Duration::from_millis(TICK_TIME)).await;
    }
}