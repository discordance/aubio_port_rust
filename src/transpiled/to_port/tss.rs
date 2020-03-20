use ::libc;
extern "C" {
    #[no_mangle]
    fn calloc(_: libc::c_ulong, _: libc::c_ulong) -> *mut libc::c_void;
    #[no_mangle]
    fn free(_: *mut libc::c_void);
    #[no_mangle]
    fn fabsf(_: libc::c_float) -> libc::c_float;
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
    /* * compute the principal argument

  This function maps the input phase to its corresponding value wrapped in the
range \f$ [-\pi, \pi] \f$.

  \param phase unwrapped phase to map to the unit circle

  \return equivalent phase wrapped to the unit circle

*/
    #[no_mangle]
    fn aubio_unwrap2pi(phase: smpl_t) -> smpl_t;
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
pub type C2RustUnnamed = libc::c_uint;
pub const AUBIO_FAIL: C2RustUnnamed = 1;
pub const AUBIO_OK: C2RustUnnamed = 0;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct cvec_t {
    pub length: uint_t,
    pub norm: *mut smpl_t,
    pub phas: *mut smpl_t,
}
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
/* default values : alpha=4, beta=3, threshold=0.25 */
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_tss_t {
    pub threshold: smpl_t,
    pub alpha: smpl_t,
    pub beta: smpl_t,
    pub parm: smpl_t,
    pub thrsfact: smpl_t,
    pub theta1: *mut fvec_t,
    pub theta2: *mut fvec_t,
    pub oft1: *mut fvec_t,
    pub oft2: *mut fvec_t,
    pub dev: *mut fvec_t,
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

  Transient / Steady-state Separation (TSS)

  This file implement a Transient / Steady-state Separation (TSS) as described
  in:

  Christopher Duxbury, Mike E. Davies, and Mark B. Sandler. Separation of
  transient information in musical audio using multiresolution analysis
  techniques. In Proceedings of the Digital Audio Effects Conference, DAFx-01,
  pages 1--5, Limerick, Ireland, 2001.

  Available at http://www.csis.ul.ie/dafx01/proceedings/papers/duxbury.pdf

  \example spectral/test-tss.c

*/
/* * Transient / Steady-state Separation object */
pub type aubio_tss_t = _aubio_tss_t;
/* * split input into transient and steady states components

  \param o tss object as returned by new_aubio_tss()
  \param input input spectral frame
  \param trans output transient components
  \param stead output steady state components

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_tss_do(mut o: *mut aubio_tss_t,
                                      mut input: *const cvec_t,
                                      mut trans: *mut cvec_t,
                                      mut stead: *mut cvec_t) {
    let mut j: uint_t = 0;
    let mut test: uint_t = 0;
    let mut nbins: uint_t = (*input).length;
    let mut alpha: smpl_t = (*o).alpha;
    let mut beta: smpl_t = (*o).beta;
    let mut parm: smpl_t = (*o).parm;
    let mut dev: *mut smpl_t = (*(*o).dev).data;
    let mut oft1: *mut smpl_t = (*(*o).oft1).data;
    let mut oft2: *mut smpl_t = (*(*o).oft2).data;
    let mut theta1: *mut smpl_t = (*(*o).theta1).data;
    let mut theta2: *mut smpl_t = (*(*o).theta2).data;
    /* second phase derivative */
    j = 0 as libc::c_int as uint_t;
    while j < nbins {
        *dev.offset(j as isize) =
            aubio_unwrap2pi((*(*input).phas.offset(j as isize) as
                                 libc::c_double -
                                 2.0f64 *
                                     *theta1.offset(j as isize) as
                                         libc::c_double +
                                 *theta2.offset(j as isize) as libc::c_double)
                                as smpl_t);
        *theta2.offset(j as isize) = *theta1.offset(j as isize);
        *theta1.offset(j as isize) = *(*input).phas.offset(j as isize);
        /* transient analysis */
        test =
            (fabsf(*dev.offset(j as isize)) > parm * *oft1.offset(j as isize))
                as libc::c_int as uint_t;
        *(*trans).norm.offset(j as isize) =
            *(*input).norm.offset(j as isize) * test as libc::c_float;
        *(*trans).phas.offset(j as isize) =
            *(*input).phas.offset(j as isize) * test as libc::c_float;
        /* steady state analysis */
        test =
            (fabsf(*dev.offset(j as isize)) < parm * *oft2.offset(j as isize))
                as libc::c_int as uint_t;
        *(*stead).norm.offset(j as isize) =
            *(*input).norm.offset(j as isize) * test as libc::c_float;
        *(*stead).phas.offset(j as isize) =
            *(*input).phas.offset(j as isize) * test as libc::c_float;
        /*increase probability for transient */
        test =
            (*(*trans).norm.offset(j as isize) as libc::c_double == 0.0f64) as
                libc::c_int as uint_t;
        *oft1.offset(j as isize) = test as smpl_t;
        test =
            (*(*trans).norm.offset(j as isize) as libc::c_double > 0.0f64) as
                libc::c_int as uint_t;
        let ref mut fresh0 = *oft1.offset(j as isize);
        *fresh0 += alpha * test as libc::c_float;
        test =
            (*oft1.offset(j as isize) as libc::c_double > 1.0f64 &&
                 *(*trans).norm.offset(j as isize) as libc::c_double > 0.0f64)
                as libc::c_int as uint_t;
        let ref mut fresh1 = *oft1.offset(j as isize);
        *fresh1 += beta * test as libc::c_float;
        /*increase probability for steady states */
        test =
            (*(*stead).norm.offset(j as isize) as libc::c_double == 0.0f64) as
                libc::c_int as uint_t;
        *oft2.offset(j as isize) = test as smpl_t;
        test =
            (*(*stead).norm.offset(j as isize) as libc::c_double > 0.0f64) as
                libc::c_int as uint_t;
        let ref mut fresh2 = *oft2.offset(j as isize);
        *fresh2 += alpha * test as libc::c_float;
        test =
            (*oft2.offset(j as isize) as libc::c_double > 1.0f64 &&
                 *(*stead).norm.offset(j as isize) as libc::c_double > 0.0f64)
                as libc::c_int as uint_t;
        let ref mut fresh3 = *oft2.offset(j as isize);
        *fresh3 += beta * test as libc::c_float;
        j = j.wrapping_add(1)
    };
}
/* * set transient / steady state separation threshold

  \param o tss object as returned by new_aubio_tss()
  \param thrs new threshold value

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_tss_set_threshold(mut o: *mut aubio_tss_t,
                                                 mut threshold: smpl_t)
 -> uint_t {
    (*o).threshold = threshold;
    (*o).parm = (*o).threshold * (*o).thrsfact;
    return AUBIO_OK as libc::c_int as uint_t;
}
/* * create tss object

  \param buf_size buffer size
  \param hop_size step size

*/
#[no_mangle]
pub unsafe extern "C" fn new_aubio_tss(mut buf_size: uint_t,
                                       mut hop_size: uint_t)
 -> *mut aubio_tss_t {
    let mut o: *mut aubio_tss_t =
        calloc(::std::mem::size_of::<aubio_tss_t>() as libc::c_ulong,
               1 as libc::c_int as libc::c_ulong) as *mut aubio_tss_t;
    let mut rsize: uint_t =
        buf_size.wrapping_div(2 as libc::c_int as
                                  libc::c_uint).wrapping_add(1 as libc::c_int
                                                                 as
                                                                 libc::c_uint);
    (*o).threshold = 0.25f64 as smpl_t;
    (*o).thrsfact =
        (3.14159265358979323846264338327950288f64 * 2.0f64 *
             hop_size as libc::c_double / rsize as libc::c_double) as smpl_t;
    (*o).alpha = 3.0f64 as smpl_t;
    (*o).beta = 4.0f64 as smpl_t;
    (*o).parm = (*o).threshold * (*o).thrsfact;
    (*o).theta1 = new_fvec(rsize);
    (*o).theta2 = new_fvec(rsize);
    (*o).oft1 = new_fvec(rsize);
    (*o).oft2 = new_fvec(rsize);
    (*o).dev = new_fvec(rsize);
    return o;
}
/* * delete tss object

  \param o tss object as returned by new_aubio_tss()

*/
#[no_mangle]
pub unsafe extern "C" fn del_aubio_tss(mut s: *mut aubio_tss_t) {
    del_fvec((*s).theta1);
    del_fvec((*s).theta2);
    del_fvec((*s).oft1);
    del_fvec((*s).oft2);
    del_fvec((*s).dev);
    free(s as *mut libc::c_void);
}
/* * set parameter a, defaults to 3

  \param o tss object as returned by new_aubio_tss()
  \param alpha new value for alpha parameter

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_tss_set_alpha(mut o: *mut aubio_tss_t,
                                             mut alpha: smpl_t) -> uint_t {
    (*o).alpha = alpha;
    return AUBIO_OK as libc::c_int as uint_t;
}
/* * set parameter b, defaults to 3

  \param o tss object as returned by new_aubio_tss()
  \param beta new value for beta parameter

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_tss_set_beta(mut o: *mut aubio_tss_t,
                                            mut beta: smpl_t) -> uint_t {
    (*o).beta = beta;
    return AUBIO_OK as libc::c_int as uint_t;
}
