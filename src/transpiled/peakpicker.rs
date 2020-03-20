use crate::transpiled::filter::aubio_filter_t;

extern "C" {
    // pub type _aubio_filter_t;
    #[no_mangle]
    fn free(_: *mut core::ffi::c_void);
    #[no_mangle]
    fn calloc(_: u64, _: u64) -> *mut core::ffi::c_void;
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
    /* * make a copy of a vector

  \param s source vector
  \param t vector to copy to

*/
    #[no_mangle]
    fn fvec_copy(s: *const fvec_t, t: *mut fvec_t);
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

  Various math functions

  \example test-mathutils.c
  \example test-mathutils-window.c

 */
    /* * compute the mean of a vector

  \param s vector to compute mean from
  \return the mean of `v`

*/
    #[no_mangle]
    fn fvec_mean(s: *mut fvec_t) -> smpl_t;
    /* * push a new element to the end of a vector, erasing the first element and
 * sliding all others

  \param in vector to push to
  \param new_elem new_element to add at the end of the vector

  In numpy words, this is equivalent to: in = np.concatenate([in, [new_elem]])[1:]

*/
    #[no_mangle]
    fn fvec_push(in_0: *mut fvec_t, new_elem: smpl_t);
    /* * returns the median of a vector

  The QuickSelect routine is based on the algorithm described in "Numerical
recipes in C", Second Edition, Cambridge University Press, 1992, Section 8.5,
ISBN 0-521-43108-5

  This implementation of the QuickSelect routine is based on Nicolas
Devillard's implementation, available at http://ndevilla.free.fr/median/median/
and in the Public Domain.

  \param v vector to get median from

  \return the median of v

*/
    #[no_mangle]
    fn fvec_median(v: *mut fvec_t) -> smpl_t;
    /* * finds exact peak index by quadratic interpolation

  See [Quadratic Interpolation of Spectral
  Peaks](https://ccrma.stanford.edu/~jos/sasp/Quadratic_Peak_Interpolation.html),
  by Julius O. Smith III

  \f$ p_{frac} = \frac{1}{2} \frac {x[p-1] - x[p+1]} {x[p-1] - 2 x[p] + x[p+1]} \in [ -.5, .5] \f$

  \param x vector to get the interpolated peak position from
  \param p index of the peak in vector `x`
  \return \f$ p + p_{frac} \f$ exact peak position of interpolated maximum or minimum

*/
    #[no_mangle]
    fn fvec_quadratic_peak_pos(x: *const fvec_t, p: uint_t) -> smpl_t;
    /* * return 1 if v[p] is a peak and positive, 0 otherwise

  This function returns 1 if a peak is found at index p in the vector v. The
peak is defined as follows:

  - v[p] is positive
  - v[p-1] < v[p]
  - v[p] > v[p+1]

  \param v input vector
  \param p position of supposed for peak

  \return 1 if a peak is found, 0 otherwise

*/
    #[no_mangle]
    fn fvec_peakpick(v: *const fvec_t, p: uint_t) -> uint_t;
    /* * filter input vector forward and backward

  \param f ::aubio_filter_t object as returned by new_aubio_filter()
  \param in ::fvec_t input vector to filter
  \param tmp memory space to use for computation

*/
    #[no_mangle]
    fn aubio_filter_do_filtfilt(f: *mut aubio_filter_t, in_0: *mut fvec_t,
                                tmp: *mut fvec_t);
    /* * delete a filter object

  \param f filter object to delete

*/
    #[no_mangle]
    fn del_aubio_filter(f: *mut aubio_filter_t);
    /* * create biquad filter with `b0`, `b1`, `b2`, `a1`, `a2` coeffs

  \param b0 forward filter coefficient
  \param b1 forward filter coefficient
  \param b2 forward filter coefficient
  \param a1 feedback filter coefficient
  \param a2 feedback filter coefficient

*/
    #[no_mangle]
    fn new_aubio_filter_biquad(b0: lsmp_t, b1: lsmp_t, b2: lsmp_t, a1: lsmp_t,
                               a2: lsmp_t) -> *mut aubio_filter_t;
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
pub type smpl_t = f32;
/* * print format for sample in single precision */
/* * long sample format (64 bits or more) */
pub type lsmp_t = f64;
/* * print format for sample in double precision */
/* * unsigned integer */
pub type uint_t = u32;
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
pub type C2RustUnnamed = u32;
pub const AUBIO_FAIL: C2RustUnnamed = 1;
pub const AUBIO_OK: C2RustUnnamed = 0;
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

  Digital filter

  This object stores a digital filter of order \f$n\f$.
  It contains the following data:
    - \f$ n*1 b_i \f$ feedforward coefficients
    - \f$ n*1 a_i \f$ feedback coefficients
    - \f$ n*c x_i \f$ input signal
    - \f$ n*c y_i \f$ output signal

  For convenience, the samplerate of the input signal is also stored in the
  object.

  Feedforward and feedback parameters can be modified using
  aubio_filter_get_feedback() and aubio_filter_get_feedforward().

  The function aubio_filter_do_outplace() computes the following output signal
  \f$ y[n] \f$ from the input signal \f$ x[n] \f$:

  \f{eqnarray*}{
     y[n] = b_0 x[n] & + & b_1 x[n-1] + b_2 x[n-2] + ... + b_P x[n-P] \\
                     & - & a_1 y[n-1] - a_2 y[n-2] - ... - a_P y[n-P] \\
  \f}

  The function aubio_filter_do() executes the same computation but modifies
  directly the input signal (in-place).

  The function aubio_filter_do_filtfilt() version runs the filter twice, first
  forward then backward, to compensate with the phase shifting of the forward
  operation.

  Some convenience functions are provided:
    - new_aubio_filter_a_weighting() and aubio_filter_set_a_weighting(),
    - new_aubio_filter_c_weighting() and aubio_filter_set_c_weighting().
    - new_aubio_filter_biquad() and aubio_filter_set_biquad().

  \example temporal/test-filter.c

*/
/* peak picking parameters, default values in brackets
 *
 *     [<----post----|--pre-->]
 *  .................|.............
 *  time->           ^now
 */
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_peakpicker_t {
    pub threshold: smpl_t,
    pub win_post: uint_t,
    pub win_pre: uint_t,
    pub thresholdfn: aubio_thresholdfn_t,
    pub pickerfn: aubio_pickerfn_t,
    pub biquad: *mut aubio_filter_t,
    pub onset_keep: *mut fvec_t,
    pub onset_proc: *mut fvec_t,
    pub onset_peek: *mut fvec_t,
    pub thresholded: *mut fvec_t,
    pub scratch: *mut fvec_t,
}
/* * function pointer to peak-picking function */
pub type aubio_pickerfn_t
    =
    Option<unsafe extern "C" fn(_: *mut fvec_t, _: uint_t) -> uint_t>;
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
/* * function pointer to thresholding function */
pub type aubio_thresholdfn_t
    =
    Option<unsafe extern "C" fn(_: *mut fvec_t) -> smpl_t>;
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

  Peak picking utilities function

  \example onset/test-peakpicker.c

*/
/* * peak-picker structure */
pub type aubio_peakpicker_t = _aubio_peakpicker_t;
/* * real time peak picking function */
/* * \bug should be used to calculate filter coefficients */
  /* cutoff: low-pass filter cutoff [0.34, 1] */
  /* smpl_t cutoff; */
/* not used anymore */
  /* time precision [512/44100  winlength/samplerate, fs/buffer_size */
  /* smpl_t tau; */
  /* alpha: normalisation exponent [9] */
  /* smpl_t alpha; */
/* * modified version for real time, moving mean adaptive threshold this method
 * is slightly more permissive than the offline one, and yelds to an increase
 * of false positives. best  */
#[no_mangle]
pub unsafe extern "C" fn aubio_peakpicker_do(mut p: *mut aubio_peakpicker_t,
                                             mut onset: *mut fvec_t,
                                             mut out: *mut fvec_t) {
    let mut onset_keep: *mut fvec_t = (*p).onset_keep;
    let mut onset_proc: *mut fvec_t = (*p).onset_proc;
    let mut onset_peek: *mut fvec_t = (*p).onset_peek;
    let mut thresholded: *mut fvec_t = (*p).thresholded;
    let mut scratch: *mut fvec_t = (*p).scratch;
    let mut mean: smpl_t = 0.0f64 as smpl_t;
    let mut median: smpl_t = 0.0f64 as smpl_t;
    let mut j: uint_t = 0 as i32 as uint_t;
    /* push new novelty to the end */
    fvec_push(onset_keep, *(*onset).data.offset(0 as i32 as isize));
    /* store a copy */
    fvec_copy(onset_keep, onset_proc);
    /* filter this copy */
    aubio_filter_do_filtfilt((*p).biquad, onset_proc, scratch);
    /* calculate mean and median for onset_proc */
    mean = fvec_mean(onset_proc);
    /* copy to scratch and compute its median */
    fvec_copy(onset_proc, scratch);
    median = (*p).thresholdfn.expect("non-null function pointer")(scratch);
    /* shift peek array */
    j = 0 as i32 as uint_t;
    while j < (3 as i32 - 1 as i32) as u32 {
        *(*onset_peek).data.offset(j as isize) =
            *(*onset_peek).data.offset(j.wrapping_add(1 as i32 as
                                                          u32) as
                                           isize);
        j = j.wrapping_add(1)
    }
    /* calculate new tresholded value */
    *(*thresholded).data.offset(0 as i32 as isize) =
        *(*onset_proc).data.offset((*p).win_post as isize) - median -
            mean * (*p).threshold;
    *(*onset_peek).data.offset(2 as i32 as isize) =
        *(*thresholded).data.offset(0 as i32 as isize);
    *(*out).data.offset(0 as i32 as isize) =
        (*p).pickerfn.expect("non-null function pointer")(onset_peek,
                                                          1 as i32 as
                                                              uint_t) as
            smpl_t;
    if *(*out).data.offset(0 as i32 as isize) != 0. {
        *(*out).data.offset(0 as i32 as isize) =
            fvec_quadratic_peak_pos(onset_peek, 1 as i32 as uint_t)
    };
}
/* * get current peak value */
/* * this method returns the current value in the pick peaking buffer
 * after smoothing
 */
#[no_mangle]
pub unsafe extern "C" fn aubio_peakpicker_get_thresholded_input(mut p:
                                                                    *mut aubio_peakpicker_t)
 -> *mut fvec_t {
    return (*p).thresholded;
}
/* * set peak picking threshold */
#[no_mangle]
pub unsafe extern "C" fn aubio_peakpicker_set_threshold(mut p:
                                                            *mut aubio_peakpicker_t,
                                                        mut threshold: smpl_t)
 -> uint_t {
    (*p).threshold = threshold;
    return AUBIO_OK as i32 as uint_t;
}
/* * get peak picking threshold */
#[no_mangle]
pub unsafe extern "C" fn aubio_peakpicker_get_threshold(mut p:
                                                            *mut aubio_peakpicker_t)
 -> smpl_t {
    return (*p).threshold;
}
/* * set peak picker thresholding function */
#[no_mangle]
pub unsafe extern "C" fn aubio_peakpicker_set_thresholdfn(mut p:
                                                              *mut aubio_peakpicker_t,
                                                          mut thresholdfn:
                                                              aubio_thresholdfn_t)
 -> uint_t {
    (*p).thresholdfn = thresholdfn;
    return AUBIO_OK as i32 as uint_t;
}
/* * get peak picker thresholding function */
#[no_mangle]
pub unsafe extern "C" fn aubio_peakpicker_get_thresholdfn(mut p:
                                                              *mut aubio_peakpicker_t)
 -> aubio_thresholdfn_t {
    return (*p).thresholdfn;
}
/* * peak-picker creation function */
#[no_mangle]
pub unsafe extern "C" fn new_aubio_peakpicker() -> *mut aubio_peakpicker_t {
    let mut t: *mut aubio_peakpicker_t =
        calloc(::std::mem::size_of::<aubio_peakpicker_t>() as u64,
               1 as i32 as u64) as
            *mut aubio_peakpicker_t; /* 0.0668; 0.33; 0.082; 0.033; */
    (*t).threshold = 0.1f64 as smpl_t; /* (fvec_mean); */
    (*t).win_post = 5 as i32 as uint_t;
    (*t).win_pre = 1 as i32 as uint_t;
    (*t).thresholdfn =
        ::std::mem::transmute::<Option<unsafe extern "C" fn(_: *mut fvec_t)
                                           -> smpl_t>,
                                aubio_thresholdfn_t>(Some(fvec_median as
                                                              unsafe extern "C" fn(_:
                                                                                       *mut fvec_t)
                                                                  -> smpl_t));
    (*t).pickerfn =
        ::std::mem::transmute::<Option<unsafe extern "C" fn(_: *const fvec_t,
                                                            _: uint_t)
                                           -> uint_t>,
                                aubio_pickerfn_t>(Some(fvec_peakpick as
                                                           unsafe extern "C" fn(_:
                                                                                    *const fvec_t,
                                                                                _:
                                                                                    uint_t)
                                                               -> uint_t));
    (*t).scratch =
        new_fvec((*t).win_post.wrapping_add((*t).win_pre).wrapping_add(1 as
                                                                           i32
                                                                           as
                                                                           u32));
    (*t).onset_keep =
        new_fvec((*t).win_post.wrapping_add((*t).win_pre).wrapping_add(1 as
                                                                           i32
                                                                           as
                                                                           u32));
    (*t).onset_proc =
        new_fvec((*t).win_post.wrapping_add((*t).win_pre).wrapping_add(1 as
                                                                           i32
                                                                           as
                                                                           u32));
    (*t).onset_peek = new_fvec(3 as i32 as uint_t);
    (*t).thresholded = new_fvec(1 as i32 as uint_t);
    /* cutoff: low-pass filter with cutoff reduced frequency at 0.34
     generated with octave butter function: [b,a] = butter(2, 0.34);
   */
    (*t).biquad =
        new_aubio_filter_biquad(0.15998789f64, 0.31997577f64, 0.15998789f64,
                                0.23484048f64, 0 as i32 as lsmp_t);
    return t;
}
/* * destroy peak picker structure */
#[no_mangle]
pub unsafe extern "C" fn del_aubio_peakpicker(mut p:
                                                  *mut aubio_peakpicker_t) {
    del_aubio_filter((*p).biquad);
    del_fvec((*p).onset_keep);
    del_fvec((*p).onset_proc);
    del_fvec((*p).onset_peek);
    del_fvec((*p).thresholded);
    del_fvec((*p).scratch);
    free(p as *mut core::ffi::c_void);
}
