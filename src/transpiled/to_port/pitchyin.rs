use ::libc;
extern "C" {
    #[no_mangle]
    fn calloc(_: libc::c_ulong, _: libc::c_ulong) -> *mut libc::c_void;
    #[no_mangle]
    fn free(_: *mut libc::c_void);
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
    /* * find the index of the min of a vector

  \param s vector to get the index from

  \return the index of the minimum element of v

*/
    #[no_mangle]
    fn fvec_min_elem(s: *mut fvec_t) -> uint_t;
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
/* This algorithm was developed by A. de Cheveigné and H. Kawahara and
 * published in:
 * 
 * de Cheveigné, A., Kawahara, H. (2002) "YIN, a fundamental frequency
 * estimator for speech and music", J. Acoust. Soc. Am. 111, 1917-1930.  
 *
 * see http://recherche.ircam.fr/equipes/pcm/pub/people/cheveign.html
 */
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_pitchyin_t {
    pub yin: *mut fvec_t,
    pub tol: smpl_t,
    pub peak_pos: uint_t,
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

  Pitch detection using the YIN algorithm

  This algorithm was developed by A. de Cheveigne and H. Kawahara and
  published in:

  De Cheveigné, A., Kawahara, H. (2002) "YIN, a fundamental frequency
  estimator for speech and music", J. Acoust. Soc. Am. 111, 1917-1930.

  see http://recherche.ircam.fr/equipes/pcm/pub/people/cheveign.html
      http://recherche.ircam.fr/equipes/pcm/cheveign/ps/2002_JASA_YIN_proof.pdf

  \example pitch/test-pitchyin.c

*/
/* * pitch detection object */
pub type aubio_pitchyin_t = _aubio_pitchyin_t;
/* * creation of the pitch detection object

  \param buf_size size of the input buffer to analyse

*/
#[no_mangle]
pub unsafe extern "C" fn new_aubio_pitchyin(mut bufsize: uint_t)
 -> *mut aubio_pitchyin_t {
    let mut o: *mut aubio_pitchyin_t =
        calloc(::std::mem::size_of::<aubio_pitchyin_t>() as libc::c_ulong,
               1 as libc::c_int as libc::c_ulong) as *mut aubio_pitchyin_t;
    (*o).yin =
        new_fvec(bufsize.wrapping_div(2 as libc::c_int as libc::c_uint));
    (*o).tol = 0.15f64 as smpl_t;
    (*o).peak_pos = 0 as libc::c_int as uint_t;
    return o;
}
/* * deletion of the pitch detection object

  \param o pitch detection object as returned by new_aubio_pitchyin()

*/
#[no_mangle]
pub unsafe extern "C" fn del_aubio_pitchyin(mut o: *mut aubio_pitchyin_t) {
    del_fvec((*o).yin);
    free(o as *mut libc::c_void);
}
/* * execute pitch detection an input buffer

  \param o pitch detection object as returned by new_aubio_pitchyin()
  \param samples_in input signal vector (length as specified at creation time)
  \param cands_out pitch period candidates, in samples

*/
/* all the above in one */
#[no_mangle]
pub unsafe extern "C" fn aubio_pitchyin_do(mut o: *mut aubio_pitchyin_t,
                                           mut input: *const fvec_t,
                                           mut out: *mut fvec_t) {
    let tol: smpl_t = (*o).tol;
    let mut yin: *mut fvec_t = (*o).yin;
    let mut input_data: *const smpl_t = (*input).data;
    let length: uint_t = (*yin).length;
    let mut yin_data: *mut smpl_t = (*yin).data;
    let mut j: uint_t = 0;
    let mut tau: uint_t = 0;
    let mut period: sint_t = 0;
    let mut tmp: smpl_t = 0.;
    let mut tmp2: smpl_t = 0.0f64 as smpl_t;
    *yin_data.offset(0 as libc::c_int as isize) = 1.0f64 as smpl_t;
    tau = 1 as libc::c_int as uint_t;
    while tau < length {
        *yin_data.offset(tau as isize) = 0.0f64 as smpl_t;
        j = 0 as libc::c_int as uint_t;
        while j < length {
            tmp =
                *input_data.offset(j as isize) -
                    *input_data.offset(j.wrapping_add(tau) as isize);
            let ref mut fresh0 = *yin_data.offset(tau as isize);
            *fresh0 += tmp * tmp;
            j = j.wrapping_add(1)
        }
        tmp2 += *yin_data.offset(tau as isize);
        if tmp2 != 0 as libc::c_int as libc::c_float {
            let ref mut fresh1 = *(*yin).data.offset(tau as isize);
            *fresh1 *= tau as libc::c_float / tmp2
        } else { *(*yin).data.offset(tau as isize) = 1.0f64 as smpl_t }
        period = tau.wrapping_sub(3 as libc::c_int as libc::c_uint) as sint_t;
        if tau > 4 as libc::c_int as libc::c_uint &&
               *yin_data.offset(period as isize) < tol &&
               *yin_data.offset(period as isize) <
                   *yin_data.offset((period + 1 as libc::c_int) as isize) {
            (*o).peak_pos = period as uint_t;
            *(*out).data.offset(0 as libc::c_int as isize) =
                fvec_quadratic_peak_pos(yin, (*o).peak_pos);
            return
        }
        tau = tau.wrapping_add(1)
    }
    (*o).peak_pos = fvec_min_elem(yin);
    *(*out).data.offset(0 as libc::c_int as isize) =
        fvec_quadratic_peak_pos(yin, (*o).peak_pos);
}
/* * get current confidence of YIN algorithm

  \param o YIN pitch detection object
  \return confidence parameter

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_pitchyin_get_confidence(mut o:
                                                           *mut aubio_pitchyin_t)
 -> smpl_t {
    return (1.0f64 -
                *(*(*o).yin).data.offset((*o).peak_pos as isize) as
                    libc::c_double) as smpl_t;
}
/* * set tolerance parameter for YIN algorithm

  \param o YIN pitch detection object
  \param tol tolerance parameter for minima selection [default 0.15]

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_pitchyin_set_tolerance(mut o:
                                                          *mut aubio_pitchyin_t,
                                                      mut tol: smpl_t)
 -> uint_t {
    (*o).tol = tol;
    return 0 as libc::c_int as uint_t;
}
/* * get tolerance parameter for YIN algorithm

  \param o YIN pitch detection object
  \return tolerance parameter for minima selection [default 0.15]

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_pitchyin_get_tolerance(mut o:
                                                          *mut aubio_pitchyin_t)
 -> smpl_t {
    return (*o).tol;
}
