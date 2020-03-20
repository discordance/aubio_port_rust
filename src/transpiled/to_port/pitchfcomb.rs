use ::libc;
extern "C" {
    pub type _aubio_fft_t;
    #[no_mangle]
    fn calloc(_: libc::c_ulong, _: libc::c_ulong) -> *mut libc::c_void;
    #[no_mangle]
    fn free(_: *mut libc::c_void);
    #[no_mangle]
    fn log10f(_: libc::c_float) -> libc::c_float;
    #[no_mangle]
    fn memmove(_: *mut libc::c_void, _: *const libc::c_void, _: libc::c_ulong)
     -> *mut libc::c_void;
    /* * fvec_t buffer creation function

  \param length the length of the buffer to create

*/
    #[no_mangle]
    fn new_fvec(length: uint_t) -> *mut fvec_t;
    /* * fvec_t buffer deletion function

  \param s buffer to delete as returned by new_fvec()

*/
    #[no_mangle]
    fn del_fvec(s: *mut fvec_t);
    #[no_mangle]
    fn new_aubio_window(window_type: *mut char_t, size: uint_t)
     -> *mut fvec_t;
    #[no_mangle]
    fn aubio_unwrap2pi(phase: smpl_t) -> smpl_t;
    /* * cvec_t buffer creation function

  This function creates a cvec_t structure holding two arrays of size
  [length/2+1], corresponding to the norm and phase values of the
  spectral frame. The length stored in the structure is the actual size of both
  arrays, not the length of the complex and symmetrical vector, specified as
  creation argument.

  \param length the length of the buffer to create

*/
    #[no_mangle]
    fn new_cvec(length: uint_t) -> *mut cvec_t;
    /* * cvec_t buffer deletion function

  \param s buffer to delete as returned by new_cvec()

*/
    #[no_mangle]
    fn del_cvec(s: *mut cvec_t);
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
}
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
pub type smpl_t = libc::c_float;
/* * print format for sample in double precision */
/* * unsigned integer */
pub type uint_t = libc::c_uint;
/* * signed integer */
pub type sint_t = libc::c_int;
/* * character */
pub type char_t = libc::c_char;
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

  Vector of real-valued data

  This file specifies the ::fvec_t buffer type, which is used throughout aubio
  to store vector of real-valued ::smpl_t.

  \example test-fvec.c

*/
/* * Buffer for real data

  Vector of real-valued data

  ::fvec_t is is the structure used to store vector of real-valued data, ::smpl_t .

  \code

  uint_t buffer_size = 1024;

  // create a vector of 512 values
  fvec_t * input = new_fvec (buffer_size);

  // set some values of the vector
  input->data[23] = 2.;
  // ..

  // compute the mean of the vector
  mean = fvec_mean(a_vector);

  // destroy the vector
  del_fvec(a_vector);

  \endcode

  See `examples/` and `tests/src` directories for more examples.

 */
#[derive(Copy, Clone)]
#[repr(C)]
pub struct fvec_t {
    pub length: uint_t,
    pub data: *mut smpl_t,
}
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
/* * FFT object

  This object computes forward and backward FFTs.

*/
pub type aubio_fft_t = _aubio_fft_t;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_pitchfcomb_t {
    pub fftSize: uint_t,
    pub stepSize: uint_t,
    pub rate: uint_t,
    pub winput: *mut fvec_t,
    pub win: *mut fvec_t,
    pub fftOut: *mut cvec_t,
    pub fftLastPhase: *mut fvec_t,
    pub fft: *mut aubio_fft_t,
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

   Pitch detection using a fast harmonic comb filter

   This pitch extraction method implements a fast harmonic comb filter to
   determine the fundamental frequency of a harmonic sound.

   This file was derived from the tuneit project, written by Mario Lang to
   detect the fundamental frequency of a sound.

   See http://delysid.org/tuneit.html

   \example pitch/test-pitchfcomb.c

*/
/* * pitch detection object */
pub type aubio_pitchfcomb_t = _aubio_pitchfcomb_t;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct aubio_fpeak_t {
    pub bin: smpl_t,
    pub db: smpl_t,
}
/* * creation of the pitch detection object

  \param buf_size size of the input buffer to analyse
  \param hop_size step size between two consecutive analysis instant

*/
#[no_mangle]
pub unsafe extern "C" fn new_aubio_pitchfcomb(mut bufsize: uint_t,
                                              mut hopsize: uint_t)
 -> *mut aubio_pitchfcomb_t {
    let mut p: *mut aubio_pitchfcomb_t =
        calloc(::std::mem::size_of::<aubio_pitchfcomb_t>() as libc::c_ulong,
               1 as libc::c_int as libc::c_ulong) as *mut aubio_pitchfcomb_t;
    (*p).fftSize = bufsize;
    (*p).stepSize = hopsize;
    (*p).fft = new_aubio_fft(bufsize);
    if (*p).fft.is_null() {
        free(p as *mut libc::c_void);
        return 0 as *mut aubio_pitchfcomb_t
    } else {
        (*p).winput = new_fvec(bufsize);
        (*p).fftOut = new_cvec(bufsize);
        (*p).fftLastPhase = new_fvec(bufsize);
        (*p).win =
            new_aubio_window(b"hanning\x00" as *const u8 as
                                 *const libc::c_char as *mut char_t, bufsize);
        return p
    };
}
/* * execute pitch detection on an input buffer

  \param p pitch detection object as returned by new_aubio_pitchfcomb
  \param input input signal window (length as specified at creation time)
  \param output pitch candidates in bins

*/
/* input must be stepsize long */
#[no_mangle]
pub unsafe extern "C" fn aubio_pitchfcomb_do(mut p: *mut aubio_pitchfcomb_t,
                                             mut input: *const fvec_t,
                                             mut output: *mut fvec_t) {
    let mut k: uint_t = 0;
    let mut l: uint_t = 0;
    let mut maxharm: uint_t = 0 as libc::c_int as uint_t;
    let mut phaseDifference: smpl_t =
        (3.14159265358979323846264338327950288f64 * 2.0f64 *
             (*p).stepSize as smpl_t as libc::c_double /
             (*p).fftSize as smpl_t as libc::c_double) as smpl_t;
    let mut peaks: [aubio_fpeak_t; 8] = [aubio_fpeak_t{bin: 0., db: 0.,}; 8];
    k = 0 as libc::c_int as uint_t;
    while k < 8 as libc::c_int as libc::c_uint {
        peaks[k as usize].db = -200.0f64 as smpl_t;
        peaks[k as usize].bin = 0.0f64 as smpl_t;
        k = k.wrapping_add(1)
    }
    k = 0 as libc::c_int as uint_t;
    while k < (*input).length {
        *(*(*p).winput).data.offset(k as isize) =
            *(*(*p).win).data.offset(k as isize) *
                *(*input).data.offset(k as isize);
        k = k.wrapping_add(1)
    }
    aubio_fft_do((*p).fft, (*p).winput, (*p).fftOut);
    k = 0 as libc::c_int as uint_t;
    while k <= (*p).fftSize.wrapping_div(2 as libc::c_int as libc::c_uint) {
        let mut magnitude: smpl_t =
            (20.0f64 *
                 log10f((2.0f64 *
                             *(*(*p).fftOut).norm.offset(k as isize) as
                                 libc::c_double /
                             (*p).fftSize as smpl_t as libc::c_double) as
                            libc::c_float) as libc::c_double) as smpl_t;
        let mut phase: smpl_t = *(*(*p).fftOut).phas.offset(k as isize);
        let mut tmp: smpl_t = 0.;
        let mut bin: smpl_t = 0.;
        /* compute phase difference */
        tmp = phase - *(*(*p).fftLastPhase).data.offset(k as isize);
        *(*(*p).fftLastPhase).data.offset(k as isize) = phase;
        /* subtract expected phase difference */
        tmp -= k as smpl_t * phaseDifference;
        /* map delta phase into +/- Pi interval */
        tmp = aubio_unwrap2pi(tmp);
        /* get deviation from bin frequency from the +/- Pi interval */
        tmp =
            (((*p).fftSize as libc::c_float / (*p).stepSize as smpl_t * tmp)
                 as libc::c_double /
                 (3.14159265358979323846264338327950288f64 * 2.0f64)) as
                smpl_t;
        /* compute the k-th partials' true bin */
        bin = k as smpl_t + tmp;
        if bin as libc::c_double > 0.0f64 &&
               magnitude > peaks[0 as libc::c_int as usize].db {
            // && magnitude < 0) {
            memmove(peaks.as_mut_ptr().offset(1 as libc::c_int as isize) as
                        *mut libc::c_void,
                    peaks.as_mut_ptr() as *const libc::c_void,
                    (::std::mem::size_of::<aubio_fpeak_t>() as
                         libc::c_ulong).wrapping_mul((8 as libc::c_int -
                                                          1 as libc::c_int) as
                                                         libc::c_ulong));
            peaks[0 as libc::c_int as usize].bin = bin;
            peaks[0 as libc::c_int as usize].db = magnitude
        }
        k = k.wrapping_add(1)
    }
    k = 0 as libc::c_int as uint_t;
    l = 1 as libc::c_int as uint_t;
    while l < 8 as libc::c_int as libc::c_uint &&
              peaks[l as usize].bin as libc::c_double > 0.0f64 {
        let mut harmonic: sint_t = 0;
        harmonic = 5 as libc::c_int;
        while harmonic > 1 as libc::c_int {
            if ((peaks[0 as libc::c_int as usize].bin / peaks[l as usize].bin)
                    as libc::c_double) < harmonic as libc::c_double + 0.02f64
                   &&
                   (peaks[0 as libc::c_int as usize].bin /
                        peaks[l as usize].bin) as libc::c_double >
                       harmonic as libc::c_double - 0.02f64 {
                if harmonic > maxharm as sint_t &&
                       peaks[0 as libc::c_int as usize].db <
                           peaks[l as usize].db /
                               2 as libc::c_int as libc::c_float {
                    maxharm = harmonic as uint_t;
                    k = l
                }
            }
            harmonic -= 1
        }
        l = l.wrapping_add(1)
    }
    *(*output).data.offset(0 as libc::c_int as isize) = peaks[k as usize].bin;
    /* quick hack to clean output a bit */
    if peaks[k as usize].bin as libc::c_double > 5000.0f64 {
        *(*output).data.offset(0 as libc::c_int as isize) = 0.0f64 as smpl_t
    };
}
/* * deletion of the pitch detection object

  \param p pitch detection object as returned by new_aubio_pitchfcomb

*/
#[no_mangle]
pub unsafe extern "C" fn del_aubio_pitchfcomb(mut p:
                                                  *mut aubio_pitchfcomb_t) {
    del_cvec((*p).fftOut);
    del_fvec((*p).fftLastPhase);
    del_fvec((*p).win);
    del_fvec((*p).winput);
    del_aubio_fft((*p).fft);
    free(p as *mut libc::c_void);
}
