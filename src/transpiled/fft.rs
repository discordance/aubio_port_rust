

extern "C" {
    #[no_mangle]
    fn atan2f(_: f32, _: f32) -> f32;
    #[no_mangle]
    fn cosf(_: f32) -> f32;
    #[no_mangle]
    fn sinf(_: f32) -> f32;
    #[no_mangle]
    fn calloc(_: u64, _: u64) -> *mut core::ffi::c_void;
    #[no_mangle]
    fn free(_: *mut core::ffi::c_void);
    #[no_mangle]
    fn aubio_log(level: sint_t, fmt: *const char_t, _: ...) -> uint_t;
    #[no_mangle]
    fn aubio_is_power_of_two(a: uint_t) -> uint_t;
    #[no_mangle]
    fn del_fvec(s: *mut fvec_t);
    #[no_mangle]
    fn new_fvec(length: uint_t) -> *mut fvec_t;
    #[no_mangle]
    fn fabsf(_: f32) -> f32;
    #[no_mangle]
    fn sqrtf(_: f32) -> f32;
    /*
  Copyright (C) 2003-2009 Paul Brossier <piem@aubio.org>

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
    // using FFTW3
    // using ACCELERATE
    // using INTEL IPP
    // using OOURA
    // let's use ooura instead
    #[no_mangle]
    fn aubio_ooura_rdft(_: i32, _: i32, _: *mut smpl_t,
                        _: *mut i32, _: *mut smpl_t);
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
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_fft_t {
    pub winsize: uint_t,
    pub fft_size: uint_t,
    pub in_0: *mut smpl_t,
    pub out: *mut smpl_t,
    pub w: *mut smpl_t,
    pub ip: *mut i32,
    pub compspec: *mut fvec_t,
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
/* * FFT object

  This object computes forward and backward FFTs.

*/
pub type aubio_fft_t = _aubio_fft_t;
/* * create new FFT computation object

  \param size length of the FFT

*/
#[no_mangle]
pub unsafe extern "C" fn new_aubio_fft(winsize: uint_t)
 -> *mut aubio_fft_t {
    let mut s: *mut aubio_fft_t =
        calloc(::std::mem::size_of::<aubio_fft_t>() as u64,
               1 as i32 as u64) as *mut aubio_fft_t;
    if (winsize as sint_t) < 2 as i32 {
        aubio_log(AUBIO_LOG_ERR as i32,
                  b"AUBIO ERROR: fft: got winsize %d, but can not be < 2\n\x00"
                      as *const u8 as *const i8, winsize);
    } else if aubio_is_power_of_two(winsize) !=
                  1 as i32 as u32 {
        aubio_log(AUBIO_LOG_ERR as i32,
                  b"AUBIO ERROR: fft: can only create with sizes power of two, requested %d, try recompiling aubio with --enable-fftw3\n\x00"
                      as *const u8 as *const i8, winsize);
    } else {
        (*s).winsize = winsize;
        (*s).fft_size =
            winsize.wrapping_div(2 as i32 as
                                     u32).wrapping_add(1 as
                                                                    i32
                                                                    as
                                                                    u32);
        (*s).compspec = new_fvec(winsize);
        (*s).in_0 =
            calloc(((*s).winsize as
                        u64).wrapping_mul(::std::mem::size_of::<smpl_t>()
                                                        as u64),
                   1 as i32 as u64) as *mut smpl_t;
        (*s).out =
            calloc(((*s).winsize as
                        u64).wrapping_mul(::std::mem::size_of::<smpl_t>()
                                                        as u64),
                   1 as i32 as u64) as *mut smpl_t;
        (*s).ip =
            calloc(((*s).fft_size as
                        u64).wrapping_mul(::std::mem::size_of::<i32>()
                                                        as u64),
                   1 as i32 as u64) as *mut i32;
        (*s).w =
            calloc(((*s).fft_size as
                        u64).wrapping_mul(::std::mem::size_of::<smpl_t>()
                                                        as u64),
                   1 as i32 as u64) as *mut smpl_t;
        *(*s).ip.offset(0 as i32 as isize) = 0 as i32;
        // using ACCELERATE
        // using Intel IPP
        // using OOURA
        /* using OOURA */
        return s
    }
    free(s as *mut core::ffi::c_void);
    return 0 as *mut aubio_fft_t;
}
/* * delete FFT object

  \param s fft object as returned by new_aubio_fft

*/
#[no_mangle]
pub unsafe extern "C" fn del_aubio_fft(s: *mut aubio_fft_t) {
    /* destroy data */
    // using FFTW3
    // using ACCELERATE
    // using Intel IPP
    // using OOURA
    free((*s).w as *mut core::ffi::c_void);
    free((*s).ip as *mut core::ffi::c_void);
    del_fvec((*s).compspec);
    free((*s).in_0 as *mut core::ffi::c_void);
    free((*s).out as *mut core::ffi::c_void);
    free(s as *mut core::ffi::c_void);
}
/* * compute forward FFT

  \param s fft object as returned by new_aubio_fft
  \param input input signal
  \param spectrum output spectrum

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_fft_do(s: *mut aubio_fft_t,
                                      input: *const fvec_t,
                                      spectrum: *mut cvec_t) {
    aubio_fft_do_complex(s, input, (*s).compspec);
    aubio_fft_get_spectrum((*s).compspec, spectrum);
}
/* * compute backward (inverse) FFT

  \param s fft object as returned by new_aubio_fft
  \param spectrum input spectrum
  \param output output signal

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_fft_rdo(s: *mut aubio_fft_t,
                                       spectrum: *const cvec_t,
                                       output: *mut fvec_t) {
    aubio_fft_get_realimag(spectrum, (*s).compspec);
    aubio_fft_rdo_complex(s, (*s).compspec, output);
}
/* * compute forward FFT

  \param s fft object as returned by new_aubio_fft
  \param input real input signal
  \param compspec complex output fft real/imag

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_fft_do_complex(s: *mut aubio_fft_t,
                                              input: *const fvec_t,
                                              compspec: *mut fvec_t) {
    let mut i: uint_t = 0;
    i = 0 as i32 as uint_t;
    while i < (*s).winsize {
        *(*s).in_0.offset(i as isize) = *(*input).data.offset(i as isize);
        i = i.wrapping_add(1)
    }
    /* HAVE_MEMCPY_HACKS */
    // using FFTW3
    // using ACCELERATE
    // using Intel IPP
    // using OOURA
    aubio_ooura_rdft((*s).winsize as i32, 1 as i32, (*s).in_0,
                     (*s).ip, (*s).w);
    *(*compspec).data.offset(0 as i32 as isize) =
        *(*s).in_0.offset(0 as i32 as isize);
    *(*compspec).data.offset((*s).winsize.wrapping_div(2 as i32 as
                                                           u32) as
                                 isize) =
        *(*s).in_0.offset(1 as i32 as isize);
    i = 1 as i32 as uint_t;
    while i < (*s).fft_size.wrapping_sub(1 as i32 as u32) {
        *(*compspec).data.offset(i as isize) =
            *(*s).in_0.offset((2 as i32 as
                                   u32).wrapping_mul(i) as isize);
        *(*compspec).data.offset((*s).winsize.wrapping_sub(i) as isize) =
            -*(*s).in_0.offset((2 as i32 as
                                    u32).wrapping_mul(i).wrapping_add(1
                                                                                   as
                                                                                   i32
                                                                                   as
                                                                                   u32)
                                   as isize);
        i = i.wrapping_add(1)
    };
    /* using OOURA */
}
/* * compute backward (inverse) FFT from real/imag

  \param s fft object as returned by new_aubio_fft
  \param compspec real/imag input fft array
  \param output real output array

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_fft_rdo_complex(s: *mut aubio_fft_t,
                                               compspec: *const fvec_t,
                                               output: *mut fvec_t) {
    let mut i: uint_t = 0;
    // using ACCELERATE
    // using Intel IPP
    // using OOURA
    let scale: smpl_t =
        (2.0f64 / (*s).winsize as f64) as smpl_t;
    *(*s).out.offset(0 as i32 as isize) =
        *(*compspec).data.offset(0 as i32 as isize);
    *(*s).out.offset(1 as i32 as isize) =
        *(*compspec).data.offset((*s).winsize.wrapping_div(2 as i32 as
                                                               u32)
                                     as isize);
    i = 1 as i32 as uint_t;
    while i < (*s).fft_size.wrapping_sub(1 as i32 as u32) {
        *(*s).out.offset((2 as i32 as u32).wrapping_mul(i) as
                             isize) = *(*compspec).data.offset(i as isize);
        *(*s).out.offset((2 as i32 as
                              u32).wrapping_mul(i).wrapping_add(1 as
                                                                             i32
                                                                             as
                                                                             u32)
                             as isize) =
            -*(*compspec).data.offset((*s).winsize.wrapping_sub(i) as isize);
        i = i.wrapping_add(1)
    }
    aubio_ooura_rdft((*s).winsize as i32, -(1 as i32),
                     (*s).out, (*s).ip, (*s).w);
    i = 0 as i32 as uint_t;
    while i < (*s).winsize {
        *(*output).data.offset(i as isize) =
            *(*s).out.offset(i as isize) * scale;
        i = i.wrapping_add(1)
    };
}
/* * convert real/imag spectrum to norm/phas spectrum

  \param compspec real/imag input fft array
  \param spectrum cvec norm/phas output array

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_fft_get_spectrum(compspec: *const fvec_t,
                                                spectrum: *mut cvec_t) {
    aubio_fft_get_phas(compspec, spectrum);
    aubio_fft_get_norm(compspec, spectrum);
}
/* * convert real/imag spectrum to norm/phas spectrum

  \param compspec real/imag input fft array
  \param spectrum cvec norm/phas output array

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_fft_get_realimag(spectrum: *const cvec_t,
                                                compspec: *mut fvec_t) {
    aubio_fft_get_imag(spectrum, compspec);
    aubio_fft_get_real(spectrum, compspec);
}
/* * compute phas spectrum from real/imag parts

  \param compspec real/imag input fft array
  \param spectrum cvec norm/phas output array

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_fft_get_phas(compspec: *const fvec_t,
                                            spectrum: *mut cvec_t) {
    let mut i: uint_t = 0;
    if *(*compspec).data.offset(0 as i32 as isize) <
           0 as i32 as f32 {
        *(*spectrum).phas.offset(0 as i32 as isize) =
            3.14159265358979323846264338327950288f64 as smpl_t
    } else {
        *(*spectrum).phas.offset(0 as i32 as isize) = 0.0f64 as smpl_t
    }
    i = 1 as i32 as uint_t;
    while i <
              (*spectrum).length.wrapping_sub(1 as i32 as
                                                  u32) {
        *(*spectrum).phas.offset(i as isize) =
            atan2f(*(*compspec).data.offset((*compspec).length.wrapping_sub(i)
                                                as isize),
                   *(*compspec).data.offset(i as isize));
        i = i.wrapping_add(1)
    }
    if *(*compspec).data.offset((*compspec).length.wrapping_div(2 as
                                                                    i32
                                                                    as
                                                                    u32)
                                    as isize) <
           0 as i32 as f32 {
        *(*spectrum).phas.offset((*spectrum).length.wrapping_sub(1 as
                                                                     i32
                                                                     as
                                                                     u32)
                                     as isize) =
            3.14159265358979323846264338327950288f64 as smpl_t
    } else {
        *(*spectrum).phas.offset((*spectrum).length.wrapping_sub(1 as
                                                                     i32
                                                                     as
                                                                     u32)
                                     as isize) = 0.0f64 as smpl_t
    };
}
/* * compute norm component from real/imag parts

  \param compspec real/imag input fft array
  \param spectrum cvec norm/phas output array

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_fft_get_norm(compspec: *const fvec_t,
                                            spectrum: *mut cvec_t) {
    let mut i: uint_t = 0 as i32 as uint_t;
    *(*spectrum).norm.offset(0 as i32 as isize) =
        fabsf(*(*compspec).data.offset(0 as i32 as isize));
    i = 1 as i32 as uint_t;
    while i <
              (*spectrum).length.wrapping_sub(1 as i32 as
                                                  u32) {
        *(*spectrum).norm.offset(i as isize) =
            sqrtf(*(*compspec).data.offset(i as isize) *
                      *(*compspec).data.offset(i as isize) +
                      *(*compspec).data.offset((*compspec).length.wrapping_sub(i)
                                                   as isize) *
                          *(*compspec).data.offset((*compspec).length.wrapping_sub(i)
                                                       as isize));
        i = i.wrapping_add(1)
    }
    *(*spectrum).norm.offset((*spectrum).length.wrapping_sub(1 as i32
                                                                 as
                                                                 u32)
                                 as isize) =
        fabsf(*(*compspec).data.offset((*compspec).length.wrapping_div(2 as
                                                                           i32
                                                                           as
                                                                           u32)
                                           as isize));
}
/* * compute imaginary part from the norm/phas cvec

  \param spectrum norm/phas input array
  \param compspec real/imag output fft array

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_fft_get_imag(spectrum: *const cvec_t,
                                            compspec: *mut fvec_t) {
    let mut i: uint_t = 0;
    i = 1 as i32 as uint_t;
    while i <
              (*compspec).length.wrapping_add(1 as i32 as
                                                  u32).wrapping_div(2
                                                                                 as
                                                                                 i32
                                                                                 as
                                                                                 u32)
          {
        *(*compspec).data.offset((*compspec).length.wrapping_sub(i) as isize)
            =
            *(*spectrum).norm.offset(i as isize) *
                sinf(*(*spectrum).phas.offset(i as isize));
        i = i.wrapping_add(1)
    };
}
/* * compute real part from norm/phas components

  \param spectrum norm/phas input array
  \param compspec real/imag output fft array

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_fft_get_real(spectrum: *const cvec_t,
                                            compspec: *mut fvec_t) {
    let mut i: uint_t = 0;
    i = 0 as i32 as uint_t;
    while i <
              (*compspec).length.wrapping_div(2 as i32 as
                                                  u32).wrapping_add(1
                                                                                 as
                                                                                 i32
                                                                                 as
                                                                                 u32)
          {
        *(*compspec).data.offset(i as isize) =
            *(*spectrum).norm.offset(i as isize) *
                cosf(*(*spectrum).phas.offset(i as isize));
        i = i.wrapping_add(1)
    };
}
