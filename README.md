# Aubio port in Rust via c2rust

Porting some Aubio library modules in pure (unsafe) Rust via __c2rust__.
The goal is to be as standalone as it can be (no libc / no static or dymamic lib)

The transpiled code is ugly unsafe Rust but it should work like the origin lib.
We try to provide an idiomatic rust wrapper around the libs.

Rust wrappers until now:

- [X] Onset detector
- [X] Tempo(BPM) detector
- [X] Phase Vocoder
- [ ] Everything else ...

## Use with cargo

```toml
[dependencies.aubio_port_rs]
git = "https://github.com/discordance/aubio_port_rust.git"
```

## Example: Onset detector

```rust
use aubio_port_rs::onset::{OnsetMode, Onset};

// settings used for onset detection
const HOP_SIZE: usize = 512;
const WIND_SIZE: usize = 2048;
const SR: usize = 44_100;

// load some 44100 signal
let mono: Vec<f32> = samples.load(...);

// chunk iterator yields HOP_SIZE samples      
let mut chunk_iter = mono.chunks(HOP_SIZE);

// onset detector
let mut onset = Onset::new(OnsetMode::SpecFlux(), WIND_SIZE, HOP_SIZE, SR).expect("Onset::new");

// some detection params
onset.set_threshold(0.3);
onset.set_silence(-30.0);
onset.set_minioi(0.02);

// hold detected positions
let mut positions: Vec<usize> = Vec::new();

// track the last detected onset in sample
let mut last_detection = 0;
loop {
    let next = chunk_iter.next();
    match next {
        Some(chunk) => {
            // break the loop because we do not have enough samples
            if chunk.len() != HOP_SIZE {
                break;
            }

            // call the onset detector aubio style
            onset.execute(&chunk);
            let detected = onset.last_onset();

            // check for some invalid bug from aubio
            if detected > mono.len() as u32 {
                continue;
            }

            // update
            if last_detection < detected {
                positions.push(detected as usize);
                last_detection = detected;
            }
        }
        None => break,
    }
}
```



