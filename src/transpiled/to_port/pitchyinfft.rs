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
    #[no_mangle]
    fn new_aubio_window(window_type: *mut char_t, size: uint_t)
     -> *mut fvec_t;
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
    /* *< length of buffer */
    /* *< data vector of length ::fvec_t.length */
    /* * fvec_t buffer creation function

  \param length the length of the buffer to create

*/
    /* * fvec_t buffer deletion function

  \param s buffer to delete as returned by new_fvec()

*/
    /* * read sample value in a buffer

  \param s vector to read from
  \param position sample position to read from

*/
    /* * write sample value in a buffer

  \param s vector to write to
  \param data value to write in s->data[position]
  \param position sample position to write to

*/
    /* * read data from a buffer

  \param s vector to read from

*/
    /* * print out fvec data

  \param s vector to print out

*/
    /* * set all elements to a given value

  \param s vector to modify
  \param val value to set elements to

*/
    /* * set all elements to zero

  \param s vector to modify

*/
    /* * set all elements to ones

  \param s vector to modify

*/
    /* * revert order of vector elements

  \param s vector to revert

*/
    /* * apply weight to vector

  If the weight vector is longer than s, only the first elements are used. If
  the weight vector is shorter than s, the last elements of s are not weighted.

  \param s vector to weight
  \param weight weighting coefficients

*/
    /* * make a copy of a vector

  \param s source vector
  \param t vector to copy to

*/
    /* * make a copy of a vector, applying weights to each element

  \param in input vector
  \param weight weights vector
  \param out output vector

*/
    #[no_mangle]
    fn fvec_weighted_copy(in_0: *const fvec_t, weight_0: *const fvec_t,
                          out: *mut fvec_t);
    #[no_mangle]
    fn del_fvec(s: *mut fvec_t);
    #[no_mangle]
    fn new_fvec(length: uint_t) -> *mut fvec_t;
    #[no_mangle]
    fn powf(_: libc::c_float, _: libc::c_float) -> libc::c_float;
    #[no_mangle]
    fn floorf(_: libc::c_float) -> libc::c_float;
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
pub type smpl_t = libc::c_float;
pub type uint_t = libc::c_uint;
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
pub type char_t = libc::c_char;
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
/* * pitch yinfft structure */
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_pitchyinfft_t {
    pub win: *mut fvec_t,
    pub winput: *mut fvec_t,
    pub sqrmag: *mut fvec_t,
    pub weight: *mut fvec_t,
    pub fftout: *mut fvec_t,
    pub fft: *mut aubio_fft_t,
    pub yinfft: *mut fvec_t,
    pub tol: smpl_t,
    pub peak_pos: uint_t,
    pub short_period: uint_t,
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

  Pitch detection using a spectral implementation of the YIN algorithm

  This algorithm was derived from the YIN algorithm. In this implementation, a
  Fourier transform is used to compute a tapered square difference function,
  which allows spectral weighting. Because the difference function is tapered,
  the selection of the period is simplified.

  Paul Brossier, [Automatic annotation of musical audio for interactive
  systems](http://aubio.org/phd/), Chapter 3, Pitch Analysis, PhD thesis,
  Centre for Digital music, Queen Mary University of London, London, UK, 2006.

  \example pitch/test-pitchyinfft.c

*/
/* * pitch detection object */
pub type aubio_pitchyinfft_t = _aubio_pitchyinfft_t;
/* * shortest period under which to check for octave error */
static mut freqs: [smpl_t; 35] =
    [0.0f64 as smpl_t, 20.0f64 as smpl_t, 25.0f64 as smpl_t,
     31.5f64 as smpl_t, 40.0f64 as smpl_t, 50.0f64 as smpl_t,
     63.0f64 as smpl_t, 80.0f64 as smpl_t, 100.0f64 as smpl_t,
     125.0f64 as smpl_t, 160.0f64 as smpl_t, 200.0f64 as smpl_t,
     250.0f64 as smpl_t, 315.0f64 as smpl_t, 400.0f64 as smpl_t,
     500.0f64 as smpl_t, 630.0f64 as smpl_t, 800.0f64 as smpl_t,
     1000.0f64 as smpl_t, 1250.0f64 as smpl_t, 1600.0f64 as smpl_t,
     2000.0f64 as smpl_t, 2500.0f64 as smpl_t, 3150.0f64 as smpl_t,
     4000.0f64 as smpl_t, 5000.0f64 as smpl_t, 6300.0f64 as smpl_t,
     8000.0f64 as smpl_t, 9000.0f64 as smpl_t, 10000.0f64 as smpl_t,
     12500.0f64 as smpl_t, 15000.0f64 as smpl_t, 20000.0f64 as smpl_t,
     25100.0f64 as smpl_t, -1.0f64 as smpl_t];
static mut weight: [smpl_t; 34] =
    [-75.8f64 as smpl_t, -70.1f64 as smpl_t, -60.8f64 as smpl_t,
     -52.1f64 as smpl_t, -44.2f64 as smpl_t, -37.5f64 as smpl_t,
     -31.3f64 as smpl_t, -25.6f64 as smpl_t, -20.9f64 as smpl_t,
     -16.5f64 as smpl_t, -12.6f64 as smpl_t, -9.60f64 as smpl_t,
     -7.00f64 as smpl_t, -4.70f64 as smpl_t, -3.00f64 as smpl_t,
     -1.80f64 as smpl_t, -0.80f64 as smpl_t, -0.20f64 as smpl_t,
     -0.00f64 as smpl_t, 0.50f64 as smpl_t, 1.60f64 as smpl_t,
     3.20f64 as smpl_t, 5.40f64 as smpl_t, 7.80f64 as smpl_t,
     8.10f64 as smpl_t, 5.30f64 as smpl_t, -2.40f64 as smpl_t,
     -11.1f64 as smpl_t, -12.8f64 as smpl_t, -12.2f64 as smpl_t,
     -7.40f64 as smpl_t, -17.8f64 as smpl_t, -17.8f64 as smpl_t,
     -17.8f64 as smpl_t];
/* * creation of the pitch detection object

  \param samplerate samplerate of the input signal
  \param buf_size size of the input buffer to analyse

*/
#[no_mangle]
pub unsafe extern "C" fn new_aubio_pitchyinfft(mut samplerate: uint_t,
                                               mut bufsize: uint_t)
 -> *mut aubio_pitchyinfft_t {
    let mut i: uint_t = 0 as libc::c_int as uint_t;
    let mut j: uint_t = 1 as libc::c_int as uint_t;
    let mut freq: smpl_t = 0 as libc::c_int as smpl_t;
    let mut a0: smpl_t = 0 as libc::c_int as smpl_t;
    let mut a1: smpl_t = 0 as libc::c_int as smpl_t;
    let mut f0: smpl_t = 0 as libc::c_int as smpl_t;
    let mut f1: smpl_t = 0 as libc::c_int as smpl_t;
    let mut p: *mut aubio_pitchyinfft_t =
        calloc(::std::mem::size_of::<aubio_pitchyinfft_t>() as libc::c_ulong,
               1 as libc::c_int as libc::c_ulong) as *mut aubio_pitchyinfft_t;
    (*p).winput = new_fvec(bufsize);
    (*p).fft = new_aubio_fft(bufsize);
    if (*p).fft.is_null() {
        if !(*p).winput.is_null() { del_fvec((*p).winput); }
        free(p as *mut libc::c_void);
        return 0 as *mut aubio_pitchyinfft_t
    } else {
        (*p).fftout = new_fvec(bufsize);
        (*p).sqrmag = new_fvec(bufsize);
        (*p).yinfft =
            new_fvec(bufsize.wrapping_div(2 as libc::c_int as
                                              libc::c_uint).wrapping_add(1 as
                                                                             libc::c_int
                                                                             as
                                                                             libc::c_uint));
        (*p).tol = 0.85f64 as smpl_t;
        (*p).peak_pos = 0 as libc::c_int as uint_t;
        (*p).win =
            new_aubio_window(b"hanningz\x00" as *const u8 as
                                 *const libc::c_char as *mut char_t, bufsize);
        (*p).weight =
            new_fvec(bufsize.wrapping_div(2 as libc::c_int as
                                              libc::c_uint).wrapping_add(1 as
                                                                             libc::c_int
                                                                             as
                                                                             libc::c_uint));
        i = 0 as libc::c_int as uint_t;
        while i < (*(*p).weight).length {
            freq = i as smpl_t / bufsize as smpl_t * samplerate as smpl_t;
            while freq > freqs[j as usize] &&
                      freqs[j as usize] > 0 as libc::c_int as libc::c_float {
                //p->weight->data[i] = SQRT(DB2LIN(p->weight->data[i]));
                //AUBIO_DBG("freq %3.5f > %3.5f \tsamplerate %d (Hz) \t"
      //    "(weight length %d, bufsize %d) %d %d\n", freq, freqs[j],
      //    samplerate, p->weight->length, bufsize, i, j);
                j =
                    (j as
                         libc::c_uint).wrapping_add(1 as libc::c_int as
                                                        libc::c_uint) as
                        uint_t as uint_t
            }
            a0 =
                weight[j.wrapping_sub(1 as libc::c_int as libc::c_uint) as
                           usize];
            f0 =
                freqs[j.wrapping_sub(1 as libc::c_int as libc::c_uint) as
                          usize];
            a1 = weight[j as usize];
            f1 = freqs[j as usize];
            if f0 == f1 {
                // just in case
                *(*(*p).weight).data.offset(i as isize) = a0
            } else if f0 == 0 as libc::c_int as libc::c_float {
                // y = ax+b
                *(*(*p).weight).data.offset(i as isize) =
                    (a1 - a0) / f1 * freq + a0
            } else {
                *(*(*p).weight).data.offset(i as isize) =
                    (((a1 - a0) / (f1 - f0) * freq) as libc::c_double +
                         (a0 as libc::c_double -
                              (a1 - a0) as libc::c_double /
                                  ((f1 / f0) as libc::c_double - 1.0f64))) as
                        smpl_t
            }
            while freq > freqs[j as usize] {
                j =
                    (j as
                         libc::c_uint).wrapping_add(1 as libc::c_int as
                                                        libc::c_uint) as
                        uint_t as uint_t
            }
            *(*(*p).weight).data.offset(i as isize) =
                powf(10.0f64 as libc::c_float,
                     *(*(*p).weight).data.offset(i as isize) * 0.05f32);
            i = i.wrapping_add(1)
        }
        //AUBIO_DBG("%f\n",p->weight->data[i]);
        // check for octave errors above 1300 Hz
        (*p).short_period =
            floorf((samplerate as libc::c_double / 1300.0f64 + 0.5f64) as
                       libc::c_float) as uint_t;
        return p
    };
}
/* * execute pitch detection on an input buffer

  \param o pitch detection object as returned by new_aubio_pitchyinfft
  \param samples_in input signal vector (length as specified at creation time)
  \param cands_out pitch period candidates, in samples

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_pitchyinfft_do(mut p: *mut aubio_pitchyinfft_t,
                                              mut input: *const fvec_t,
                                              mut output: *mut fvec_t) {
    let mut tau: uint_t = 0;
    let mut l: uint_t = 0;
    let mut length: uint_t = (*(*p).fftout).length;
    let mut halfperiod: uint_t = 0;
    let mut fftout: *mut fvec_t = (*p).fftout;
    let mut yin: *mut fvec_t = (*p).yinfft;
    let mut tmp: smpl_t = 0.0f64 as smpl_t;
    let mut sum: smpl_t = 0.0f64 as smpl_t;
    // window the input
    fvec_weighted_copy(input, (*p).win, (*p).winput);
    // get the real / imag parts of its fft
    aubio_fft_do_complex((*p).fft, (*p).winput, fftout);
    // get the squared magnitude spectrum, applying some weight
    *(*(*p).sqrmag).data.offset(0 as libc::c_int as isize) =
        *(*fftout).data.offset(0 as libc::c_int as isize) *
            *(*fftout).data.offset(0 as libc::c_int as isize);
    let ref mut fresh0 =
        *(*(*p).sqrmag).data.offset(0 as libc::c_int as isize);
    *fresh0 *= *(*(*p).weight).data.offset(0 as libc::c_int as isize);
    l = 1 as libc::c_int as uint_t;
    while l < length.wrapping_div(2 as libc::c_int as libc::c_uint) {
        *(*(*p).sqrmag).data.offset(l as isize) =
            *(*fftout).data.offset(l as isize) *
                *(*fftout).data.offset(l as isize) +
                *(*fftout).data.offset(length.wrapping_sub(l) as isize) *
                    *(*fftout).data.offset(length.wrapping_sub(l) as isize);
        let ref mut fresh1 = *(*(*p).sqrmag).data.offset(l as isize);
        *fresh1 *= *(*(*p).weight).data.offset(l as isize);
        *(*(*p).sqrmag).data.offset(length.wrapping_sub(l) as isize) =
            *(*(*p).sqrmag).data.offset(l as isize);
        l = l.wrapping_add(1)
    }
    *(*(*p).sqrmag).data.offset(length.wrapping_div(2 as libc::c_int as
                                                        libc::c_uint) as
                                    isize) =
        *(*fftout).data.offset(length.wrapping_div(2 as libc::c_int as
                                                       libc::c_uint) as isize)
            *
            *(*fftout).data.offset(length.wrapping_div(2 as libc::c_int as
                                                           libc::c_uint) as
                                       isize);
    let ref mut fresh2 =
        *(*(*p).sqrmag).data.offset(length.wrapping_div(2 as libc::c_int as
                                                            libc::c_uint) as
                                        isize);
    *fresh2 *=
        *(*(*p).weight).data.offset(length.wrapping_div(2 as libc::c_int as
                                                            libc::c_uint) as
                                        isize);
    // get sum of weighted squared mags
    l = 0 as libc::c_int as uint_t;
    while l <
              length.wrapping_div(2 as libc::c_int as
                                      libc::c_uint).wrapping_add(1 as
                                                                     libc::c_int
                                                                     as
                                                                     libc::c_uint)
          {
        sum += *(*(*p).sqrmag).data.offset(l as isize);
        l = l.wrapping_add(1)
    }
    sum = (sum as libc::c_double * 2.0f64) as smpl_t;
    // get the real / imag parts of the fft of the squared magnitude
    aubio_fft_do_complex((*p).fft, (*p).sqrmag, fftout);
    *(*yin).data.offset(0 as libc::c_int as isize) = 1.0f64 as smpl_t;
    tau = 1 as libc::c_int as uint_t;
    while tau < (*yin).length {
        // compute the square differences
        *(*yin).data.offset(tau as isize) =
            sum - *(*fftout).data.offset(tau as isize);
        // and the cumulative mean normalized difference function
        tmp += *(*yin).data.offset(tau as isize);
        if tmp != 0 as libc::c_int as libc::c_float {
            let ref mut fresh3 = *(*yin).data.offset(tau as isize);
            *fresh3 *= tau as libc::c_float / tmp
        } else { *(*yin).data.offset(tau as isize) = 1.0f64 as smpl_t }
        tau = tau.wrapping_add(1)
    }
    // find best candidates
    tau = fvec_min_elem(yin);
    if *(*yin).data.offset(tau as isize) < (*p).tol {
        // no interpolation, directly return the period as an integer
    //output->data[0] = tau;
    //return;
        // 3 point quadratic interpolation
    //return fvec_quadratic_peak_pos (yin,tau,1);
    /* additional check for (unlikely) octave doubling in higher frequencies */
        if tau > (*p).short_period {
            *(*output).data.offset(0 as libc::c_int as isize) =
                fvec_quadratic_peak_pos(yin, tau)
        } else {
            /* should compare the minimum value of each interpolated peaks */
            halfperiod =
                floorf((tau.wrapping_div(2 as libc::c_int as libc::c_uint) as
                            libc::c_double + 0.5f64) as libc::c_float) as
                    uint_t;
            if *(*yin).data.offset(halfperiod as isize) < (*p).tol {
                (*p).peak_pos = halfperiod
            } else { (*p).peak_pos = tau }
            *(*output).data.offset(0 as libc::c_int as isize) =
                fvec_quadratic_peak_pos(yin, (*p).peak_pos)
        }
    } else {
        (*p).peak_pos = 0 as libc::c_int as uint_t;
        *(*output).data.offset(0 as libc::c_int as isize) = 0.0f64 as smpl_t
    };
}
/* * deletion of the pitch detection object

  \param o pitch detection object as returned by new_aubio_pitchyinfft()

*/
#[no_mangle]
pub unsafe extern "C" fn del_aubio_pitchyinfft(mut p:
                                                   *mut aubio_pitchyinfft_t) {
    del_fvec((*p).win);
    del_aubio_fft((*p).fft);
    del_fvec((*p).yinfft);
    del_fvec((*p).sqrmag);
    del_fvec((*p).fftout);
    del_fvec((*p).winput);
    del_fvec((*p).weight);
    free(p as *mut libc::c_void);
}
/* * get current confidence of YIN algorithm

  \param o YIN pitch detection object
  \return confidence parameter

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_pitchyinfft_get_confidence(mut o:
                                                              *mut aubio_pitchyinfft_t)
 -> smpl_t {
    return (1.0f64 -
                *(*(*o).yinfft).data.offset((*o).peak_pos as isize) as
                    libc::c_double) as smpl_t;
}
/* * set tolerance parameter for YIN algorithm

  \param o YIN pitch detection object
  \param tol tolerance parameter for minima selection [default 0.15]

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_pitchyinfft_set_tolerance(mut p:
                                                             *mut aubio_pitchyinfft_t,
                                                         mut tol: smpl_t)
 -> uint_t {
    (*p).tol = tol;
    return 0 as libc::c_int as uint_t;
}
/* * get tolerance parameter for YIN algorithm

  \param o YIN pitch detection object

  \return tolerance parameter for minima selection [default 0.15]

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_pitchyinfft_get_tolerance(mut p:
                                                             *mut aubio_pitchyinfft_t)
 -> smpl_t {
    return (*p).tol;
}
