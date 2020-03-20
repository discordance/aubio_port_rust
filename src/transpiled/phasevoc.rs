use crate::transpiled::fft::aubio_fft_t;


use crate::transpiled::log::aubio_log;

extern "C" {
    #[no_mangle]
    fn calloc(_: u64, _: u64) -> *mut core::ffi::c_void;
    #[no_mangle]
    fn free(_: *mut core::ffi::c_void);
    /* file interface */
    /* strings */
    /* Error reporting */
    /* Logging */
    /* * internal logging function, defined in utils/log.c */
    #[no_mangle]
    fn fvec_ishift(v: *mut fvec_t);
    #[no_mangle]
    fn fvec_shift(v: *mut fvec_t);
    #[no_mangle]
    fn fvec_set_window(window: *mut fvec_t, window_type: *mut char_t)
     -> uint_t;
    #[no_mangle]
    fn new_aubio_window(window_type: *mut char_t, size: uint_t)
     -> *mut fvec_t;
    #[no_mangle]
    fn fvec_weight(s: *mut fvec_t, weight: *const fvec_t);
    #[no_mangle]
    fn del_fvec(s: *mut fvec_t);
    #[no_mangle]
    fn new_fvec(length: uint_t) -> *mut fvec_t;
    /* * create new FFT computation object

  \param size length of the FFT

*/
    #[no_mangle]
    fn new_aubio_fft(size: uint_t) -> *mut aubio_fft_t;
    /* * delete FFT object

  \param s fft object as returned by new_aubio_fft

*/
    #[no_mangle]
    fn del_aubio_fft(s: *mut aubio_fft_t);
    /* * compute forward FFT

  \param s fft object as returned by new_aubio_fft
  \param input input signal
  \param spectrum output spectrum

*/
    #[no_mangle]
    fn aubio_fft_do(s: *mut aubio_fft_t, input: *const fvec_t,
                    spectrum: *mut cvec_t);
    /* * compute backward (inverse) FFT

  \param s fft object as returned by new_aubio_fft
  \param spectrum input spectrum
  \param output output signal

*/
    #[no_mangle]
    fn aubio_fft_rdo(s: *mut aubio_fft_t, spectrum: *const cvec_t,
                     output: *mut fvec_t);
}
pub type smpl_t = f32;
pub type uint_t = u32;
pub type sint_t = i32;
/*
  Copyright (C) 2003-2015 Paul Brossier <piem@aubio.org>

  This file is part of aubio.

  aubio is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  aubio is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with aubio.  If not, see <http://www.gnu.org/licenses/>.

*/
/* * \file

  Definition of data types used in aubio

*/
/* * defined to 1 if aubio is compiled in double precision */
/* * short sample format (32 or 64 bits) */
/* * print format for sample in single precision */
/* * long sample format (64 bits or more) */
/* * print format for sample in double precision */
/* * unsigned integer */
/* * signed integer */
/* * character */
pub type char_t = i8;

#[derive(Copy, Clone)]
#[repr(C)]
pub struct fvec_t {
    pub length: uint_t,
    pub data: *mut smpl_t,
}
pub type aubio_log_level = u32;
pub const AUBIO_LOG_LAST_LEVEL: aubio_log_level = 5;
pub const AUBIO_LOG_WRN: aubio_log_level = 4;
pub const AUBIO_LOG_DBG: aubio_log_level = 3;
pub const AUBIO_LOG_MSG: aubio_log_level = 2;
pub const AUBIO_LOG_INF: aubio_log_level = 1;
pub const AUBIO_LOG_ERR: aubio_log_level = 0;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct cvec_t {
    pub length: uint_t,
    pub norm: *mut smpl_t,
    pub phas: *mut smpl_t,
}
/*
  Copyright (C) 2003-2013 Paul Brossier <piem@aubio.org>

  This file is part of aubio.

  aubio is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  aubio is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with aubio.  If not, see <http://www.gnu.org/licenses/>.

*/
/* * \file

  Fast Fourier Transform

  Depending on how aubio was compiled, FFT are computed using one of:
    - [Ooura](http://www.kurims.kyoto-u.ac.jp/~ooura/fft.html)
    - [FFTW3](http://www.fftw.org)
    - [vDSP](https://developer.apple.com/library/mac/#documentation/Accelerate/Reference/vDSPRef/Reference/reference.html)

  \example spectral/test-fft.c

*/
/*
  Copyright (C) 2003-2014 Paul Brossier <piem@aubio.org>

  This file is part of aubio.

  aubio is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  aubio is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with aubio.  If not, see <http://www.gnu.org/licenses/>.

*/
/* * phasevocoder internal object */
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_pvoc_t {
    pub win_s: uint_t,
    pub hop_s: uint_t,
    pub fft: *mut aubio_fft_t,
    pub data: *mut fvec_t,
    pub dataold: *mut fvec_t,
    pub synth: *mut fvec_t,
    pub synthold: *mut fvec_t,
    pub w: *mut fvec_t,
    pub start: uint_t,
    pub end: uint_t,
    pub scale: smpl_t,
    pub end_datasize: uint_t,
    pub hop_datasize: uint_t,
}
/*
  Copyright (C) 2003-2013 Paul Brossier <piem@aubio.org>

  This file is part of aubio.

  aubio is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  aubio is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with aubio.  If not, see <http://www.gnu.org/licenses/>.

*/
/* * \file

  Phase vocoder object

  This object implements a phase vocoder. The spectral frames are computed
  using a HanningZ window and a swapped version of the signal to simplify the
  phase relationships across frames. The window sizes and overlap are specified
  at creation time.

  \example spectral/test-phasevoc.c

*/
/* * phasevocoder object */
pub type aubio_pvoc_t = _aubio_pvoc_t;
/* * compute spectral frame

  This function accepts an input vector of size [hop_s]. The
  analysis buffer is rotated and filled with the new data. After windowing of
  this signal window, the Fourier transform is computed and returned in
  fftgrain as two vectors, magnitude and phase.

  \param pv phase vocoder object as returned by new_aubio_pvoc
  \param in new input signal (hop_s long)
  \param fftgrain output spectral frame

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_pvoc_do(pv: *mut aubio_pvoc_t,
                                       datanew: *const fvec_t,
                                       fftgrain: *mut cvec_t) {
    /* slide  */
    aubio_pvoc_swapbuffers(pv, datanew);
    /* windowing */
    fvec_weight((*pv).data, (*pv).w);
    /* shift */
    fvec_shift((*pv).data);
    /* calculate fft */
    aubio_fft_do((*pv).fft, (*pv).data, fftgrain);
}
/* * compute signal from spectral frame

  This function takes an input spectral frame fftgrain of size
  [buf_s] and computes its inverse Fourier transform. Overlap-add
  synthesis is then computed using the previously synthetised frames, and the
  output stored in out.

  \param pv phase vocoder object as returned by new_aubio_pvoc
  \param fftgrain input spectral frame
  \param out output signal (hop_s long)

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_pvoc_rdo(pv: *mut aubio_pvoc_t,
                                        fftgrain: *mut cvec_t,
                                        synthnew: *mut fvec_t) {
    /* calculate rfft */
    aubio_fft_rdo((*pv).fft, fftgrain, (*pv).synth);
    /* unshift */
    fvec_ishift((*pv).synth);
    /* windowing */
  // if overlap = 50%, do not apply window (identity)
    if (*pv).hop_s.wrapping_mul(2) <
           (*pv).win_s {
        fvec_weight((*pv).synth, (*pv).w);
    }
    /* additive synthesis */
    aubio_pvoc_addsynth(pv, synthnew);
}
/* * create phase vocoder object

  \param win_s size of analysis buffer (and length the FFT transform)
  \param hop_s step size between two consecutive analysis

*/
#[no_mangle]
pub unsafe extern "C" fn new_aubio_pvoc(win_s: uint_t, hop_s: uint_t)
 -> *mut aubio_pvoc_t {
    let mut pv: *mut aubio_pvoc_t =
        calloc(::std::mem::size_of::<aubio_pvoc_t>() as u64,
               1) as *mut aubio_pvoc_t;
    /* if (win_s < 2*hop_s) {
    AUBIO_WRN("Hop size bigger than half the window size!\n");
  } */
    if (hop_s as sint_t) < 1 {
        aubio_log(AUBIO_LOG_ERR as i32,
                  b"AUBIO ERROR: pvoc: got hop_size %d, but can not be < 1\n\x00"
                      as *const u8 as *const i8, hop_s);
    } else if (win_s as sint_t) < 2 as i32 {
        aubio_log(AUBIO_LOG_ERR as i32,
                  b"AUBIO ERROR: pvoc: got buffer_size %d, but can not be < 2\n\x00"
                      as *const u8 as *const i8, win_s);
    } else if win_s < hop_s {
        aubio_log(AUBIO_LOG_ERR as i32,
                  b"AUBIO ERROR: pvoc: hop size (%d) is larger than win size (%d)\n\x00"
                      as *const u8 as *const i8, hop_s, win_s);
    } else {
        (*pv).fft = new_aubio_fft(win_s);
        if !(*pv).fft.is_null() {
            /* remember old */
            (*pv).data = new_fvec(win_s);
            (*pv).synth = new_fvec(win_s);
            /* new input output */
            if win_s > hop_s {
                (*pv).dataold = new_fvec(win_s.wrapping_sub(hop_s));
                (*pv).synthold = new_fvec(win_s.wrapping_sub(hop_s))
            } else {
                (*pv).dataold = new_fvec(1 as i32 as uint_t);
                (*pv).synthold = new_fvec(1 as i32 as uint_t)
            }
            (*pv).w =
                new_aubio_window(b"hanningz\x00" as *const u8 as
                                     *const i8 as *mut char_t,
                                 win_s);
            (*pv).hop_s = hop_s;
            (*pv).win_s = win_s;
            /* more than 50% overlap, overlap anyway */
            if win_s < (2 as i32 as u32).wrapping_mul(hop_s)
               {
                (*pv).start = 0 as i32 as uint_t
            } else {
                /* less than 50% overlap, reset latest grain trail */
                (*pv).start = win_s.wrapping_sub(hop_s).wrapping_sub(hop_s)
            }
            if win_s > hop_s {
                (*pv).end = win_s.wrapping_sub(hop_s)
            } else { (*pv).end = 0 as i32 as uint_t }
            (*pv).end_datasize =
                ((*pv).end as
                     u64).wrapping_mul(::std::mem::size_of::<smpl_t>()
                                                     as u64) as
                    uint_t;
            (*pv).hop_datasize =
                ((*pv).hop_s as
                     u64).wrapping_mul(::std::mem::size_of::<smpl_t>()
                                                     as u64) as
                    uint_t;
            // for reconstruction with 75% overlap
            if win_s == hop_s.wrapping_mul(4 as i32 as u32) {
                (*pv).scale = (2.0f64 / 3.0f64) as smpl_t
            } else if win_s ==
                          hop_s.wrapping_mul(8 as i32 as u32)
             {
                (*pv).scale = (1.0f64 / 3.0f64) as smpl_t
            } else if win_s ==
                          hop_s.wrapping_mul(2 as i32 as u32)
             {
                (*pv).scale = 1.0f64 as smpl_t
            } else { (*pv).scale = 0.5f64 as smpl_t }
            return pv
        }
    }
    free(pv as *mut core::ffi::c_void);
    return 0 as *mut aubio_pvoc_t;
}
/* * set window type

  \param pv phase vocoder to set the window type
  \param window_type a string representing a window

  \return 0 if successful, non-zero otherwise

 */
#[no_mangle]
pub unsafe extern "C" fn aubio_pvoc_set_window(pv: *mut aubio_pvoc_t,
                                               window: *const char_t)
 -> uint_t {
    return fvec_set_window((*pv).w, window as *mut char_t);
}
/* * delete phase vocoder object

  \param pv phase vocoder object as returned by new_aubio_pvoc

*/
#[no_mangle]
pub unsafe extern "C" fn del_aubio_pvoc(pv: *mut aubio_pvoc_t) {
    del_fvec((*pv).data);
    del_fvec((*pv).synth);
    del_fvec((*pv).dataold);
    del_fvec((*pv).synthold);
    del_fvec((*pv).w);
    del_aubio_fft((*pv).fft);
    free(pv as *mut core::ffi::c_void);
}
/* * size of memory to hop_s */
/* * returns data and dataold slided by hop_s */
unsafe extern "C" fn aubio_pvoc_swapbuffers(pv: *mut aubio_pvoc_t,
                                            new: *const fvec_t) {
    /* some convenience pointers */
    let data: *mut smpl_t = (*(*pv).data).data;
    let dataold: *mut smpl_t = (*(*pv).dataold).data;
    let datanew: *mut smpl_t = (*new).data;
    let mut i: uint_t = 0;
    i = 0 as i32 as uint_t;
    while i < (*pv).end {
        *data.offset(i as isize) = *dataold.offset(i as isize);
        i = i.wrapping_add(1)
    }
    i = 0 as i32 as uint_t;
    while i < (*pv).hop_s {
        *data.offset((*pv).end.wrapping_add(i) as isize) =
            *datanew.offset(i as isize);
        i = i.wrapping_add(1)
    }
    i = 0 as i32 as uint_t;
    while i < (*pv).end {
        *dataold.offset(i as isize) =
            *data.offset(i.wrapping_add((*pv).hop_s) as isize);
        i = i.wrapping_add(1)
    };
}
/* * do additive synthesis from 'old' and 'cur' */
unsafe extern "C" fn aubio_pvoc_addsynth(pv: *mut aubio_pvoc_t,
                                         synth_new: *mut fvec_t) {
    let mut i: uint_t = 0;
    /* some convenience pointers */
    let synth: *mut smpl_t = (*(*pv).synth).data;
    let synthold: *mut smpl_t = (*(*pv).synthold).data;
    let synthnew: *mut smpl_t = (*synth_new).data;
    /* put new result in synthnew */
    i = 0 as i32 as uint_t;
    while i < (*pv).hop_s {
        *synthnew.offset(i as isize) =
            *synth.offset(i as isize) * (*pv).scale;
        i = i.wrapping_add(1)
    }
    /* no overlap, nothing else to do */
    if (*pv).end == 0 as i32 as u32 { return }
    /* add new synth to old one */
    i = 0 as i32 as uint_t;
    while i < (*pv).hop_s {
        let ref mut fresh0 = *synthnew.offset(i as isize);
        *fresh0 += *synthold.offset(i as isize);
        i = i.wrapping_add(1)
    }
    /* shift synthold */
    i = 0 as i32 as uint_t;
    while i < (*pv).start {
        *synthold.offset(i as isize) =
            *synthold.offset(i.wrapping_add((*pv).hop_s) as isize);
        i = i.wrapping_add(1)
    }
    /* erase last frame in synthold */
    i = (*pv).start;
    while i < (*pv).end {
        *synthold.offset(i as isize) = 0.0f64 as smpl_t;
        i = i.wrapping_add(1)
    }
    /* additive synth */
    i = 0 as i32 as uint_t;
    while i < (*pv).end {
        let ref mut fresh1 = *synthold.offset(i as isize);
        *fresh1 +=
            *synth.offset(i.wrapping_add((*pv).hop_s) as isize) * (*pv).scale;
        i = i.wrapping_add(1)
    };
}
/* * get window size

  \param pv phase vocoder to get the window size from

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_pvoc_get_win(pv: *mut aubio_pvoc_t)
 -> uint_t {
    return (*pv).win_s;
}
/* * get hop size

  \param pv phase vocoder to get the hop size from

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_pvoc_get_hop(pv: *mut aubio_pvoc_t)
 -> uint_t {
    return (*pv).hop_s;
}
