use crate::transpiled::onset::*;

use std::ptr;

pub enum OnsetMode {
    Default(),
    OldDefault(),
    SpecFlux(),
    SpecDiff(),
    Energy(),
    Hfc(),
    Complex(),
    ComplexDomain(),
    Phase(),
    WPhase(),
    Mkl(),
    Kl(),
}

impl OnsetMode {
    fn value(&self) -> *const i8 {
        match self {
            OnsetMode::Default() => return b"default\0" as *const u8 as *const i8,
            OnsetMode::OldDefault() => return b"old_default\0" as *const u8 as *const i8,
            OnsetMode::SpecFlux() => return b"specflux\0" as *const u8 as *const i8,
            OnsetMode::SpecDiff() => return b"specdiff\0" as *const u8 as *const i8,
            OnsetMode::Energy() => return b"energy\0" as *const u8 as *const i8,
            OnsetMode::Hfc() => return b"hfc\0" as *const u8 as *const i8,
            OnsetMode::Complex() => return b"complex\0" as *const u8 as *const i8,
            OnsetMode::ComplexDomain() => return b"complex_domain\0" as *const u8 as *const i8,
            OnsetMode::Phase() => return b"phase\0" as *const u8 as *const i8,
            OnsetMode::WPhase() => return b"wphase\0" as *const u8 as *const i8,
            OnsetMode::Mkl() => return b"mkl\0" as *const u8 as *const i8,
            OnsetMode::Kl() => return b"kl\0" as *const u8 as *const i8,
        }
    }
}

#[derive(Debug)]
pub struct Onset {
    ptr: *mut aubio_onset_t,
    hop_size: usize,
}

unsafe impl Send for Onset {}

// all non-Sync methods take &mut self:
unsafe impl Sync for Onset {}

impl Onset {
    pub fn new(
        onset_mode: OnsetMode,
        window_size: usize,
        hop_size: usize,
        sample_rate: usize,
    ) -> Result<Self, ()> {
        let ptr = unsafe {
            new_aubio_onset(
                onset_mode.value(),
                window_size as u32,
                hop_size as u32,
                sample_rate as u32,
            )
        };

        if ptr == ptr::null_mut() {
            Err(())
        } else {
            Ok(Onset { ptr, hop_size })
        }
    }

    pub fn execute(&mut self, input_buffer: &[f32]) {
        assert!(input_buffer.len() == self.hop_size);

        let mut position = vec![0f32; 2];
        let mut tempo_fvec = fvec_t {
            data: position.as_mut_ptr(),
            length: position.len() as u32,
        };

        let input_fvec = fvec_t {
            data: input_buffer.as_ptr() as *mut f32,
            length: input_buffer.len() as u32,
        };

        unsafe {
            aubio_onset_do(self.ptr, &input_fvec, &mut tempo_fvec);
        }
    }

    pub fn last_onset(&self) -> u32 {
        unsafe { aubio_onset_get_last(self.ptr) }
    }

    pub fn set_threshold(&self, threshold: f32) {
        unsafe {
            aubio_onset_set_threshold(self.ptr, threshold);
        }
    }

    pub fn set_silence(&self, silence: f32) {
        unsafe {
            aubio_onset_set_silence(self.ptr, silence);
        }
    }

    pub fn set_minioi(&mut self, minioi: f32) {
        unsafe {
            aubio_onset_set_minioi_s(self.ptr, minioi);
        }
    }
}

impl Drop for Onset {
    fn drop(&mut self) {
        unsafe { del_aubio_onset(self.ptr) }
    }
}

#[cfg(test)]
mod tests {
    use crate::onset::*;
    use sample::Sample;

    #[test]
    fn onset_test() {
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

        let mut onset: Onset = Onset::new(OnsetMode::SpecFlux(), 2048, 512, 44_100).unwrap();
        onset.set_threshold(0.3);
        onset.set_silence(-30.0);
        onset.set_minioi(0.02);

        // detected positions
        let mut positions: Vec<usize> = Vec::new();

        let mut latest_detection = 0;

        let mut chunk_iter = mono.chunks(512);

        loop {
            let next = chunk_iter.next();
            match next {
                Some(chunk) => {
                    // break the fft
                    if chunk.len() != 512 {
                        break;
                    }
                    onset.execute(&chunk);
                    let detected = onset.last_onset();

                    // check for some invalid, bug from aubio
                    if detected > mono.len() as u32 {
                        continue;
                    }

                    if latest_detection < detected {
                        positions.push(detected as usize);
                        latest_detection = detected;
                    }
                }
                None => break,
            }
        }

        assert!(positions == [3471, 20774, 30267, 41963, 52475, 62995, 73000, 84355, 105419, 114929, 126880, 147888, 157294, 169092, 190032, 199805, 211449, 221678, 232413, 242196, 253875, 274879, 284700, 296137, 306544, 317278, 327821, 331141, 335209])
    }
}
