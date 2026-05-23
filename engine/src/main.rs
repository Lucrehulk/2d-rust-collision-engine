use collision_engine::engine::engine::Room;
use collision_engine::engine::config::config_data::{CELL_SIZE, TICK_TIME};
use std::{time::Duration, sync::Arc};
use tokio::{net::TcpListener, sync::{broadcast, mpsc, oneshot}, time::sleep};
use tokio_tungstenite::{accept_async, tungstenite::Message};
use futures_util::{StreamExt, SinkExt};

const ROOM_SIZE: f32 = 1024.0;

enum ServerCommand {
    Join { reply: oneshot::Sender<usize> },
    Leave { id: usize },
    Move { id: usize, dir: u8 },
    Resize { id: usize, delta: f32 },
}

#[tokio::main]
async fn main() {
    let mut world = Room::init(ROOM_SIZE);
    let listener = TcpListener::bind("127.0.0.1:8080").await.unwrap();
    let (tx, _) = broadcast::channel::<Vec<u8>>(16);
    let tx = Arc::new(tx);
    let (cmd_tx, mut cmd_rx) = mpsc::unbounded_channel::<ServerCommand>();

    let tx_clone = tx.clone();
    let cmd_tx_clone = cmd_tx.clone();

    tokio::spawn(async move {
        while let Ok((stream, _)) = listener.accept().await {
            let tx = tx_clone.clone();
            let cmd_tx = cmd_tx_clone.clone();
            tokio::spawn(async move {
                let ws_stream = accept_async(stream).await.unwrap();
                let (mut write, mut read) = ws_stream.split();
                
                let (reply_tx, reply_rx) = oneshot::channel();
                cmd_tx.send(ServerCommand::Join { reply: reply_tx }).unwrap();
                let my_id = reply_rx.await.unwrap();

                let mut init = vec![0];
                init.extend(ROOM_SIZE.to_le_bytes());
                init.extend(((ROOM_SIZE / CELL_SIZE) as u32).to_le_bytes());
                init.extend((my_id as u32).to_le_bytes());
                write.send(Message::Binary(init)).await.unwrap();

                let mut rx = tx.subscribe();
                
                tokio::spawn(async move {
                    while let Ok(packet) = rx.recv().await {
                        if write.send(Message::Binary(packet)).await.is_err() { break };
                    }
                });

                while let Some(Ok(Message::Binary(bin))) = read.next().await {
                    if bin.is_empty() { continue };
                    match bin[0] {
                        0 if bin.len() >= 2 => { 
                            let _ = cmd_tx.send(ServerCommand::Move { id: my_id, dir: bin[1] }); 
                        }
                        1 if bin.len() >= 5 => { 
                            let delta = f32::from_le_bytes(bin[1..5].try_into().unwrap());
                            let _ = cmd_tx.send(ServerCommand::Resize { id: my_id, delta });
                        }
                        _ => {}
                    }
                }
                let _ = cmd_tx.send(ServerCommand::Leave { id: my_id });
            });
        }
    });

    loop {
        while let Ok(cmd) = cmd_rx.try_recv() {
            match cmd {
                ServerCommand::Join { reply } => {
                    let id = world.create_entity(ROOM_SIZE / 2.0, ROOM_SIZE / 2.0, 10.0, 0.9, 0.0, 0.0, 5.0, 0.8, 15.0, 1, 1);
                    let _ = reply.send(id);
                }
                ServerCommand::Leave { id } => world.remove_entity(id),
                ServerCommand::Move { id, dir } => {
                    world.create_entity_movement_from_cardinal_direction(id, dir);
                }
                ServerCommand::Resize { id, delta } => {
                    if let Some(entity) = world.entities.get(id) {
                        if !entity.replace {
                            let new_rad = f32::clamp(entity.radius + delta, 5.0, 50.0);
                            world.resize_entity(id, new_rad, false);
                        }
                    }
                }
            }
        }

        world.update();

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