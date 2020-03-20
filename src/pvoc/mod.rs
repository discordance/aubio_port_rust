use crate::transpiled::phasevoc::*;
use std::ptr;

#[derive(Debug)]
pub struct Pvoc {
    ptr: *mut aubio_pvoc_t,
    hop_size: usize,
}

unsafe impl Send for Pvoc {}

unsafe impl Sync for Pvoc {}

impl Pvoc {
    pub fn new(window_size: usize, hop_size: usize) -> Result<Self, ()> {
        let ptr = unsafe { new_aubio_pvoc(window_size as u32, hop_size as u32) };
        //
        if ptr == ptr::null_mut() {
            Err(())
        } else {
            Ok(Pvoc { ptr, hop_size })
        }
    }

    // from signal
    pub fn from_signal(
        &mut self,
        input_buffer: &[f32],
        mut norm: &mut [f32],
        mut phas: &mut [f32],
    ) {
        assert!(input_buffer.len() == self.hop_size);

        // convert input buffer
        let input_fvec = fvec_t {
            data: input_buffer.as_ptr() as *mut f32,
            length: input_buffer.len() as u32,
        };

        // create complex output
        // let mut fftgrain = super::cvec_mut(&mut norm, &mut phas);
        let mut fftgrain = cvec_t {
            length: norm.len() as u32,
            norm: norm.as_mut_ptr(),
            phas: phas.as_mut_ptr(),
        };

        unsafe {
            aubio_pvoc_do(self.ptr, &input_fvec, &mut fftgrain);
        }
    }

    pub fn to_signal(&mut self, norm: &[f32], phas: &[f32], output_buffer: &mut [f32]) {
        // convert output buffer
        let mut output_fvec = fvec_t {
            data: output_buffer.as_mut_ptr(),
            length: output_buffer.len() as u32,
        };
        // create complex input
        let mut fftgrain = cvec_t {
            length: norm.len() as u32,
            norm: norm.as_ptr() as *mut f32,
            phas: phas.as_ptr() as *mut f32,
        };

        unsafe {
            aubio_pvoc_rdo(self.ptr, &mut fftgrain, &mut output_fvec);
        }
    }
}

impl Drop for Pvoc {
    fn drop(&mut self) {
        unsafe { del_aubio_pvoc(self.ptr) }
    }
}

#[cfg(test)]
mod tests {
    use crate::pvoc::*;

    #[test]
    fn pvoc_test() {
        let in_sig = vec![1.0f32, 0.5f32, -0.5f32, -1.0f32, 1.0f32, 0.5f32, -0.5f32, -1.0f32];

        let mut norms = vec![0f32; 8];
        let mut phas = vec![0f32; 8];

        let mut rec_sig = vec![0f32; 8];

        let mut pvoc: Pvoc = Pvoc::new(16, 8).unwrap();

        pvoc.from_signal(&in_sig[..], &mut norms[..], &mut phas[..]);

        pvoc.to_signal(&norms[..], &phas[..], &mut rec_sig);

        assert!(rec_sig == [-0.13088146, 0.10349357, -0.07076013, 0.03766468, -0.009245545, -0.010170773, 0.017628148, -0.011991352]);
    }
}
