use crate::transpiled::tempo::*;
use std::ptr;

#[derive(Debug)]
pub struct Tempo {
    ptr: *mut aubio_tempo_t,
    hop_size: usize,
}

unsafe impl Send for Tempo {}

// all non-Sync methods take &mut self:
unsafe impl Sync for Tempo {}

impl Tempo {
    pub fn new(window_size: usize, hop_size: usize, sample_rate: usize) -> Result<Self, ()> {
        const DEFAULT: *const i8 = b"default\0" as *const u8 as *const i8;

        let ptr = unsafe {
            new_aubio_tempo(
                DEFAULT,
                window_size as u32,
                hop_size as u32,
                sample_rate as u32,
            )
        };

        if ptr == ptr::null_mut() {
            Err(())
        } else {
            Ok(Tempo { ptr, hop_size })
        }
    }

    /// input_buffer length must equal hop_size!
    pub fn execute(&mut self, input_buffer: &[f32]) {
        assert!(input_buffer.len() == self.hop_size);

        let mut tempo = vec![0f32; 2];
        let mut tempo_fvec = fvec_t {
            data: tempo.as_mut_ptr(),
            length: tempo.len() as u32,
        };

        let input_fvec = fvec_t {
            data: input_buffer.as_ptr() as *mut f32,
            length: input_buffer.len() as u32,
        }; 

        unsafe {
            aubio_tempo_do(self.ptr, &input_fvec, &mut tempo_fvec);
        }
    }

    pub fn bpm(&self) -> Option<f32> {
        let bpm = unsafe { aubio_tempo_get_bpm(self.ptr) };

        if bpm == 0.0 {
            None
        } else {
            Some(bpm as f32)
        }
    }

    pub fn last_beat_ms(&self) -> f32 {
        unsafe { aubio_tempo_get_last_ms(self.ptr) }
    }
}

impl Drop for Tempo {
    fn drop(&mut self) {
        unsafe { del_aubio_tempo(self.ptr) }
    }
}


#[cfg(test)]
mod tests {
    use crate::tempo::*;
    use sample::Sample;

    #[test]
    fn tempo_test() {
        // load a test sig
        let reader = hound::WavReader::open("test_files/tech_s_16.wav").unwrap();

        let samples: Vec<f32> = reader
            .into_samples::<i16>()
            .filter_map(Result::ok)
            .map(i16::to_sample::<f32>)
            .collect();

        // mono version
        let mono: Vec<f32> = samples
            .iter()
            .step_by(2)
            .zip(samples.iter().step_by(2).skip(1))
            .map(|(l, r)| (l + r) / 2.0)
            .collect();

        let mut tempo : Tempo = Tempo::new(1024, 256, 44_100).unwrap();

        let mut chunk_iter = mono.chunks(256);

        loop {
            let next = chunk_iter.next();
            match next {
                Some(chunk) => {
                    // break the fft
                    if chunk.len() != 256 {
                        break;
                    }
                    tempo.execute(&chunk);
                }
                None => break,
            }
        }
    
        // read analysed
        let analysed_tempo = tempo.bpm().expect("Should have analysed a tempo").floor();

        assert!(analysed_tempo == 125.0);
    }
}

