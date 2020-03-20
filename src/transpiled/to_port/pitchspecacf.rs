use ::libc;
extern "C" {
    pub type _aubio_fft_t;
    #[no_mangle]
    fn calloc(_: libc::c_ulong, _: libc::c_ulong) -> *mut libc::c_void;
    #[no_mangle]
    fn free(_: *mut libc::c_void);
    #[no_mangle]
    fn fvec_quadratic_peak_pos(x: *const fvec_t, p: uint_t) -> smpl_t;
    #[no_mangle]
    fn fvec_min_elem(s: *mut fvec_t) -> uint_t;
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
  \param input real input signal
  \param compspec complex output fft real/imag

*/
    #[no_mangle]
    fn aubio_fft_do_complex(s: *mut aubio_fft_t, input: *const fvec_t,
                            compspec: *mut fvec_t);
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
/*
  Copyright (C) 2013 Paul Brossier <piem@aubio.org>

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
/* * pitch specacf structure */
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_pitchspecacf_t {
    pub win: *mut fvec_t,
    pub winput: *mut fvec_t,
    pub fft: *mut aubio_fft_t,
    pub fftout: *mut fvec_t,
    pub sqrmag: *mut fvec_t,
    pub acf: *mut fvec_t,
    pub tol: smpl_t,
    pub confidence: smpl_t,
}
/*
  Copyright (C) 2013 Paul Brossier <piem@aubio.org>

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

  Pitch detection using spectral auto correlation

  This algorithm implements pitch detection by computing the autocorrelation
  function as the cosine transform of the square spectral magnitudes.

  Anssi Klapuri. Qualitative and quantitative aspects in the design of
  periodicity esti- mation algorithms. In Proceedings of the European Signal
  Processing Conference (EUSIPCO), 2000.

  Paul Brossier, [Automatic annotation of musical audio for interactive
  systems](http://aubio.org/phd/), Chapter 3, Pitch Analysis, Autocorrelation,
  pp. 75-77, PhD thesis, Centre for Digital music, Queen Mary University of
  London, London, UK, 2006.

  \example pitch/test-pitchspecacf.c

*/
/* * pitch detection object */
pub type aubio_pitchspecacf_t = _aubio_pitchspecacf_t;
/* * creation of the pitch detection object

  \param buf_size size of the input buffer to analyse

*/
/* *< confidence */
#[no_mangle]
pub unsafe extern "C" fn new_aubio_pitchspecacf(mut bufsize: uint_t)
 -> *mut aubio_pitchspecacf_t {
    let mut p: *mut aubio_pitchspecacf_t =
        calloc(::std::mem::size_of::<aubio_pitchspecacf_t>() as libc::c_ulong,
               1 as libc::c_int as libc::c_ulong) as
            *mut aubio_pitchspecacf_t;
    (*p).fft = new_aubio_fft(bufsize);
    if (*p).fft.is_null() {
        free(p as *mut libc::c_void);
        return 0 as *mut aubio_pitchspecacf_t
    } else {
        (*p).win =
            new_aubio_window(b"hanningz\x00" as *const u8 as
                                 *const libc::c_char as *mut char_t, bufsize);
        (*p).winput = new_fvec(bufsize);
        (*p).fftout = new_fvec(bufsize);
        (*p).sqrmag = new_fvec(bufsize);
        (*p).acf =
            new_fvec(bufsize.wrapping_div(2 as libc::c_int as
                                              libc::c_uint).wrapping_add(1 as
                                                                             libc::c_int
                                                                             as
                                                                             libc::c_uint));
        (*p).tol = 1.0f64 as smpl_t;
        (*p).confidence = 0.0f64 as smpl_t;
        return p
    };
}
/* * execute pitch detection on an input buffer

  \param o pitch detection object as returned by new_aubio_pitchspecacf
  \param samples_in input signal vector (length as specified at creation time)
  \param cands_out pitch period candidates, in samples

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_pitchspecacf_do(mut p:
                                                   *mut aubio_pitchspecacf_t,
                                               mut input: *const fvec_t,
                                               mut output: *mut fvec_t) {
    let mut l: uint_t = 0;
    let mut tau: uint_t = 0;
    let mut fftout: *mut fvec_t = (*p).fftout;
    // window the input
    l = 0 as libc::c_int as uint_t;
    while l < (*input).length {
        *(*(*p).winput).data.offset(l as isize) =
            *(*(*p).win).data.offset(l as isize) *
                *(*input).data.offset(l as isize);
        l = l.wrapping_add(1)
    }
    // get the real / imag parts of its fft
    aubio_fft_do_complex((*p).fft, (*p).winput, fftout);
    l = 0 as libc::c_int as uint_t;
    while l <
              (*input).length.wrapping_div(2 as libc::c_int as
                                               libc::c_uint).wrapping_add(1 as
                                                                              libc::c_int
                                                                              as
                                                                              libc::c_uint)
          {
        *(*(*p).sqrmag).data.offset(l as isize) =
            *(*fftout).data.offset(l as isize) *
                *(*fftout).data.offset(l as isize);
        l = l.wrapping_add(1)
    }
    // get the real / imag parts of the fft of the squared magnitude
    aubio_fft_do_complex((*p).fft, (*p).sqrmag, fftout);
    // copy real part to acf
    l = 0 as libc::c_int as uint_t;
    while l <
              (*fftout).length.wrapping_div(2 as libc::c_int as
                                                libc::c_uint).wrapping_add(1
                                                                               as
                                                                               libc::c_int
                                                                               as
                                                                               libc::c_uint)
          {
        *(*(*p).acf).data.offset(l as isize) =
            *(*fftout).data.offset(l as isize);
        l = l.wrapping_add(1)
    }
    // get the minimum
    tau = fvec_min_elem((*p).acf);
    // get the interpolated minimum
    *(*output).data.offset(0 as libc::c_int as isize) =
        (fvec_quadratic_peak_pos((*p).acf, tau) as libc::c_double * 2.0f64) as
            smpl_t;
}
/* * deletion of the pitch detection object

  \param o pitch detection object as returned by new_aubio_pitchspecacf()

*/
#[no_mangle]
pub unsafe extern "C" fn del_aubio_pitchspecacf(mut p:
                                                    *mut aubio_pitchspecacf_t) {
    del_fvec((*p).win);
    del_fvec((*p).winput);
    del_aubio_fft((*p).fft);
    del_fvec((*p).sqrmag);
    del_fvec((*p).fftout);
    del_fvec((*p).acf);
    free(p as *mut libc::c_void);
}
/* * get currenct confidence for `specacf` pitch detection object

  \param o pitch detection object
  \return confidence parameter

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_pitchspecacf_get_confidence(mut o:
                                                               *const aubio_pitchspecacf_t)
 -> smpl_t {
    // no confidence for now
    return (*o).confidence;
}
/* * set tolerance parameter for `specacf` pitch detection object

  \param o pitch detection object
  \param tol tolerance parameter for minima selection [default 1.]

  \return `1` on error, `0` on success

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_pitchspecacf_set_tolerance(mut p:
                                                              *mut aubio_pitchspecacf_t,
                                                          mut tol: smpl_t)
 -> uint_t {
    (*p).tol = tol;
    return 0 as libc::c_int as uint_t;
}
/* * get tolerance parameter for `specacf` pitch detection object

  \param o pitch detection object

  \return tolerance parameter for minima selection [default 1.]

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_pitchspecacf_get_tolerance(mut p:
                                                              *const aubio_pitchspecacf_t)
 -> smpl_t {
    return (*p).tol;
}
