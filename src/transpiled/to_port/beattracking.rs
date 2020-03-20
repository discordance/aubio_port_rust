use ::libc;
extern "C" {
    #[no_mangle]
    fn calloc(_: libc::c_ulong, _: libc::c_ulong) -> *mut libc::c_void;
    #[no_mangle]
    fn free(_: *mut libc::c_void);
    #[no_mangle]
    fn expf(_: libc::c_float) -> libc::c_float;
    #[no_mangle]
    fn logf(_: libc::c_float) -> libc::c_float;
    #[no_mangle]
    fn fabsf(_: libc::c_float) -> libc::c_float;
    #[no_mangle]
    fn floorf(_: libc::c_float) -> libc::c_float;
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
    /* * set all elements to zero

  \param s vector to modify

*/
    #[no_mangle]
    fn fvec_zeros(s: *mut fvec_t);
    /* * set all elements to ones

  \param s vector to modify

*/
    #[no_mangle]
    fn fvec_ones(s: *mut fvec_t);
    /* * revert order of vector elements

  \param s vector to revert

*/
    #[no_mangle]
    fn fvec_rev(s: *mut fvec_t);
    /* * apply weight to vector

  If the weight vector is longer than s, only the first elements are used. If
  the weight vector is shorter than s, the last elements of s are not weighted.

  \param s vector to weight
  \param weight weighting coefficients

*/
    #[no_mangle]
    fn fvec_weight(s: *mut fvec_t, weight: *const fvec_t);
    /* * make a copy of a vector

  \param s source vector
  \param t vector to copy to

*/
    #[no_mangle]
    fn fvec_copy(s: *const fvec_t, t: *mut fvec_t);
    /* * find the index of the max of a vector

  \param s vector to get the index from

  \return the index of the maximum element of v

*/
    #[no_mangle]
    fn fvec_max_elem(s: *mut fvec_t) -> uint_t;
    /* * compute the sum of all elements of a vector

  \param v vector to compute the sum of

  \return the sum of v

*/
    #[no_mangle]
    fn fvec_sum(v: *mut fvec_t) -> smpl_t;
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
    /* * finds magnitude of peak by quadratic interpolation

  See [Quadratic Interpolation of Spectral
  Peaks](https://ccrma.stanford.edu/~jos/sasp/Quadratic_Peak_Interpolation.html),
  by Julius O. Smith III

  \param x vector to get the magnitude of the interpolated peak position from
  \param p index of the peak in vector `x`
  \return magnitude of interpolated peak

*/
    #[no_mangle]
    fn fvec_quadratic_peak_mag(x: *mut fvec_t, p: smpl_t) -> smpl_t;
    /* * compute normalised autocorrelation function

  \param input vector to compute autocorrelation from
  \param output vector to store autocorrelation function to

*/
    #[no_mangle]
    fn aubio_autocorr(input: *const fvec_t, output: *mut fvec_t);
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
#[derive(Copy, Clone)]
#[repr(C)]
pub struct _aubio_beattracking_t {
    pub hop_size: uint_t,
    pub samplerate: uint_t,
    pub rwv: *mut fvec_t,
    pub dfwv: *mut fvec_t,
    pub gwv: *mut fvec_t,
    pub phwv: *mut fvec_t,
    pub dfrev: *mut fvec_t,
    pub acf: *mut fvec_t,
    pub acfout: *mut fvec_t,
    pub phout: *mut fvec_t,
    pub timesig: uint_t,
    pub step: uint_t,
    pub rayparam: uint_t,
    pub lastbeat: smpl_t,
    pub counter: sint_t,
    pub flagstep: uint_t,
    pub g_var: smpl_t,
    pub gp: smpl_t,
    pub bp: smpl_t,
    pub rp: smpl_t,
    pub rp1: smpl_t,
    pub rp2: smpl_t,
}
/*
  Copyright (C) 2003-2015 Matthew Davies and Paul Brossier <piem@aubio.org>

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

  Beat tracking using a context dependant model

  This file implements the causal beat tracking algorithm designed by Matthew
  Davies and described in the following articles:

  Matthew E. P. Davies and Mark D. Plumbley. Causal tempo tracking of audio.
  In Proceedings of the International Symposium on Music Information Retrieval
  (ISMIR), pages 164­169, Barcelona, Spain, 2004.

  Matthew E. P. Davies, Paul Brossier, and Mark D. Plumbley. Beat tracking
  towards automatic musical accompaniment. In Proceedings of the Audio
  Engineering Society 118th Convention, Barcelona, Spain, May 2005.

  \example tempo/test-beattracking.c

*/
/* * beat tracking object */
pub type aubio_beattracking_t = _aubio_beattracking_t;
/* * create beat tracking object

  \param winlen length of the onset detection window
  \param hop_size number of onset detection samples [512]
  \param samplerate samplerate of the input signal

*/
#[no_mangle]
pub unsafe extern "C" fn new_aubio_beattracking(mut winlen: uint_t,
                                                mut hop_size: uint_t,
                                                mut samplerate: uint_t)
 -> *mut aubio_beattracking_t {
    let mut p: *mut aubio_beattracking_t =
        calloc(::std::mem::size_of::<aubio_beattracking_t>() as libc::c_ulong,
               1 as libc::c_int as libc::c_ulong) as
            *mut aubio_beattracking_t;
    let mut i: uint_t = 0 as libc::c_int as uint_t;
    /* default value for rayleigh weighting - sets preferred tempo to 120bpm */
    let mut rayparam: smpl_t =
        (60.0f64 * samplerate as libc::c_double / 120.0f64 /
             hop_size as libc::c_double) as smpl_t;
    let mut dfwvnorm: smpl_t =
        expf(logf(2.0f64 as libc::c_float) / rayparam *
                 winlen.wrapping_add(2 as libc::c_int as libc::c_uint) as
                     libc::c_float);
    /* length over which beat period is found [128] */
    let mut laglen: uint_t =
        winlen.wrapping_div(4 as libc::c_int as libc::c_uint);
    /* step increment - both in detection function samples -i.e. 11.6ms or
   * 1 onset frame [128] */
    let mut step: uint_t =
        winlen.wrapping_div(4 as libc::c_int as
                                libc::c_uint); /* 1.5 seconds */
    (*p).hop_size = hop_size; // constthresh empirically derived!
    (*p).samplerate = samplerate;
    (*p).lastbeat = 0 as libc::c_int as smpl_t;
    (*p).counter = 0 as libc::c_int;
    (*p).flagstep = 0 as libc::c_int as uint_t;
    (*p).g_var = 3.901f64 as smpl_t;
    (*p).rp = 1 as libc::c_int as smpl_t;
    (*p).gp = 0 as libc::c_int as smpl_t;
    (*p).rayparam = rayparam as uint_t;
    (*p).step = step;
    (*p).rwv = new_fvec(laglen);
    (*p).gwv = new_fvec(laglen);
    (*p).dfwv = new_fvec(winlen);
    (*p).dfrev = new_fvec(winlen);
    (*p).acf = new_fvec(winlen);
    (*p).acfout = new_fvec(laglen);
    (*p).phwv =
        new_fvec((2 as libc::c_int as libc::c_uint).wrapping_mul(laglen));
    (*p).phout = new_fvec(winlen);
    (*p).timesig = 0 as libc::c_int as uint_t;
    /* exponential weighting, dfwv = 0.5 when i =  43 */
    i = 0 as libc::c_int as uint_t;
    while i < winlen {
        *(*(*p).dfwv).data.offset(i as isize) =
            expf(logf(2.0f64 as libc::c_float) / rayparam *
                     i.wrapping_add(1 as libc::c_int as libc::c_uint) as
                         libc::c_float) / dfwvnorm;
        i = i.wrapping_add(1)
    }
    i = 0 as libc::c_int as uint_t;
    while i < laglen {
        *(*(*p).rwv).data.offset(i as isize) =
            (i as libc::c_double + 1.0f64) as smpl_t / (rayparam * rayparam) *
                expf((-((i as libc::c_double + 1.0f64) as smpl_t *
                            (i as libc::c_double + 1.0f64) as smpl_t) as
                          libc::c_double /
                          (2.0f64 * (rayparam * rayparam) as libc::c_double))
                         as libc::c_float);
        i = i.wrapping_add(1)
    }
    return p;
}
/* * delete beat tracking object

  \param p beat tracking object

*/
#[no_mangle]
pub unsafe extern "C" fn del_aubio_beattracking(mut p:
                                                    *mut aubio_beattracking_t) {
    del_fvec((*p).rwv);
    del_fvec((*p).gwv);
    del_fvec((*p).dfwv);
    del_fvec((*p).dfrev);
    del_fvec((*p).acf);
    del_fvec((*p).acfout);
    del_fvec((*p).phwv);
    del_fvec((*p).phout);
    free(p as *mut libc::c_void);
}
/* * track the beat

  \param bt beat tracking object
  \param dfframes current input detection function frame, smoothed by
  adaptive median threshold.
  \param out stored detected beat locations

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_beattracking_do(mut bt:
                                                   *mut aubio_beattracking_t,
                                               mut dfframe: *const fvec_t,
                                               mut output: *mut fvec_t) {
    let mut i: uint_t = 0;
    let mut k: uint_t = 0;
    let mut step: uint_t = (*bt).step;
    let mut laglen: uint_t = (*(*bt).rwv).length;
    let mut winlen: uint_t = (*(*bt).dfwv).length;
    let mut maxindex: uint_t = 0 as libc::c_int as uint_t;
    //number of harmonics in shift invariant comb filterbank
    let mut numelem: uint_t =
        4 as libc::c_int as uint_t; // beat alignment (step - lastbeat)
    let mut phase: smpl_t = 0.; // beat position
    let mut beat: smpl_t = 0.; // beat period
    let mut bp: smpl_t = 0.; // used to build shift invariant comb filterbank
    let mut a: uint_t = 0; // number of elements used to find beat phase
    let mut b: uint_t = 0;
    let mut kmax: uint_t = 0;
    /* copy dfframe, apply detection function weighting, and revert */
    fvec_copy(dfframe, (*bt).dfrev);
    fvec_weight((*bt).dfrev, (*bt).dfwv);
    fvec_rev((*bt).dfrev);
    /* compute autocorrelation function */
    aubio_autocorr(dfframe, (*bt).acf);
    /* if timesig is unknown, use metrically unbiased version of filterbank */
    if (*bt).timesig == 0 {
        numelem = 4 as libc::c_int as uint_t
    } else { numelem = (*bt).timesig }
    /* first and last output values are left intentionally as zero */
    fvec_zeros((*bt).acfout);
    /* compute shift invariant comb filterbank */
    i = 1 as libc::c_int as uint_t;
    while i < laglen.wrapping_sub(1 as libc::c_int as libc::c_uint) {
        a = 1 as libc::c_int as uint_t;
        while a <= numelem {
            b = 1 as libc::c_int as uint_t;
            while b < (2 as libc::c_int as libc::c_uint).wrapping_mul(a) {
                let ref mut fresh0 = *(*(*bt).acfout).data.offset(i as isize);
                *fresh0 =
                    (*fresh0 as libc::c_double +
                         *(*(*bt).acf).data.offset(i.wrapping_mul(a).wrapping_add(b).wrapping_sub(1
                                                                                                      as
                                                                                                      libc::c_int
                                                                                                      as
                                                                                                      libc::c_uint)
                                                       as isize) as
                             libc::c_double * 1.0f64 /
                             (2.0f64 * a as libc::c_double - 1.0f64)) as
                        smpl_t;
                b = b.wrapping_add(1)
            }
            a = a.wrapping_add(1)
        }
        i = i.wrapping_add(1)
    }
    /* apply Rayleigh weight */
    fvec_weight((*bt).acfout, (*bt).rwv);
    /* find non-zero Rayleigh period */
    maxindex = fvec_max_elem((*bt).acfout);
    if maxindex > 0 as libc::c_int as libc::c_uint &&
           maxindex <
               (*(*bt).acfout).length.wrapping_sub(1 as libc::c_int as
                                                       libc::c_uint) {
        (*bt).rp = fvec_quadratic_peak_pos((*bt).acfout, maxindex)
    } else { (*bt).rp = (*bt).rayparam as smpl_t }
    /* activate biased filterbank */
    aubio_beattracking_checkstate(bt);
    // debug metronome mode
    bp = (*bt).bp;
    /* end of biased filterbank */
    if bp == 0 as libc::c_int as libc::c_float { fvec_zeros(output); return }
    /* deliberate integer operation, could be set to 3 max eventually */
    kmax = floorf(winlen as libc::c_float / bp) as uint_t;
    /* initialize output */
    fvec_zeros((*bt).phout);
    i = 0 as libc::c_int as uint_t;
    while (i as libc::c_float) < bp {
        k = 0 as libc::c_int as uint_t;
        while k < kmax {
            let ref mut fresh1 = *(*(*bt).phout).data.offset(i as isize);
            *fresh1 +=
                *(*(*bt).dfrev).data.offset(i.wrapping_add(floorf(((bp *
                                                                        k as
                                                                            libc::c_float)
                                                                       as
                                                                       libc::c_double
                                                                       +
                                                                       0.5f64)
                                                                      as
                                                                      libc::c_float)
                                                               as uint_t) as
                                                isize);
            k = k.wrapping_add(1)
        }
        i = i.wrapping_add(1)
    }
    fvec_weight((*bt).phout, (*bt).phwv);
    /* find Rayleigh period */
    maxindex = fvec_max_elem((*bt).phout);
    if maxindex >= winlen.wrapping_sub(1 as libc::c_int as libc::c_uint) {
        /* AUBIO_BEAT_WARNINGS */
        phase = step as libc::c_float - (*bt).lastbeat
    } else { phase = fvec_quadratic_peak_pos((*bt).phout, maxindex) }
    /* take back one frame delay */
    phase = (phase as libc::c_double + 1.0f64) as smpl_t;
    // debug metronome mode
    /* reset output */
    fvec_zeros(output);
    i = 1 as libc::c_int as uint_t;
    beat = bp - phase;
    // AUBIO_DBG ("bp: %f, phase: %f, lastbeat: %f, step: %d, winlen: %d\n",
  //    bp, phase, bt->lastbeat, step, winlen);
    /* the next beat will be earlier than 60% of the tempo period
    skip this one */
    if ((step as libc::c_float - (*bt).lastbeat - phase) as libc::c_double) <
           -0.40f64 * bp as libc::c_double {
        /* AUBIO_BEAT_WARNINGS */
        beat += bp
    }
    /* start counting the beats */
    while beat + bp < 0 as libc::c_int as libc::c_float { beat += bp }
    if beat >= 0 as libc::c_int as libc::c_float {
        //AUBIO_DBG ("beat: %d, %f, %f\n", i, bp, beat);
        *(*output).data.offset(i as isize) = beat;
        i = i.wrapping_add(1)
    }
    while beat + bp <= step as libc::c_float {
        beat += bp;
        //AUBIO_DBG ("beat: %d, %f, %f\n", i, bp, beat);
        *(*output).data.offset(i as isize) = beat;
        i = i.wrapping_add(1)
    }
    (*bt).lastbeat = beat;
    /* store the number of beats in this frame as the first element */
    *(*output).data.offset(0 as libc::c_int as isize) = i as smpl_t;
}
/*
  Copyright (C) 2005-2009 Matthew Davies and Paul Brossier <piem@aubio.org>

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
/* * define to 1 to print out tracking difficulties */
#[no_mangle]
pub unsafe extern "C" fn fvec_gettimesig(mut acf: *mut fvec_t,
                                         mut acflen: uint_t, mut gp: uint_t)
 -> uint_t {
    let mut k: sint_t = 0 as libc::c_int;
    let mut three_energy: smpl_t = 0.0f64 as smpl_t;
    let mut four_energy: smpl_t = 0.0f64 as smpl_t;
    if gp < 2 as libc::c_int as libc::c_uint {
        return 4 as libc::c_int as uint_t
    }
    if acflen >
           (6 as libc::c_int as
                libc::c_uint).wrapping_mul(gp).wrapping_add(2 as libc::c_int
                                                                as
                                                                libc::c_uint)
       {
        k = -(2 as libc::c_int);
        while k < 2 as libc::c_int {
            three_energy +=
                *(*acf).data.offset((3 as libc::c_int as
                                         libc::c_uint).wrapping_mul(gp).wrapping_add(k
                                                                                         as
                                                                                         libc::c_uint)
                                        as isize);
            four_energy +=
                *(*acf).data.offset((4 as libc::c_int as
                                         libc::c_uint).wrapping_mul(gp).wrapping_add(k
                                                                                         as
                                                                                         libc::c_uint)
                                        as isize);
            k += 1
        }
    } else {
        /*Expanded to be more accurate in time sig estimation */
        k = -(2 as libc::c_int);
        while k < 2 as libc::c_int {
            three_energy +=
                *(*acf).data.offset((3 as libc::c_int as
                                         libc::c_uint).wrapping_mul(gp).wrapping_add(k
                                                                                         as
                                                                                         libc::c_uint)
                                        as isize) +
                    *(*acf).data.offset((6 as libc::c_int as
                                             libc::c_uint).wrapping_mul(gp).wrapping_add(k
                                                                                             as
                                                                                             libc::c_uint)
                                            as isize);
            four_energy +=
                *(*acf).data.offset((4 as libc::c_int as
                                         libc::c_uint).wrapping_mul(gp).wrapping_add(k
                                                                                         as
                                                                                         libc::c_uint)
                                        as isize) +
                    *(*acf).data.offset((2 as libc::c_int as
                                             libc::c_uint).wrapping_mul(gp).wrapping_add(k
                                                                                             as
                                                                                             libc::c_uint)
                                            as isize);
            k += 1
        }
    }
    return if three_energy > four_energy {
               3 as libc::c_int
           } else { 4 as libc::c_int } as uint_t;
}
#[no_mangle]
pub unsafe extern "C" fn aubio_beattracking_checkstate(mut bt:
                                                           *mut aubio_beattracking_t) {
    let mut i: uint_t = 0;
    let mut j: uint_t = 0;
    let mut a: uint_t = 0;
    let mut b: uint_t = 0;
    let mut flagconst: uint_t = 0 as libc::c_int as uint_t;
    let mut counter: sint_t = (*bt).counter;
    let mut flagstep: uint_t = (*bt).flagstep;
    let mut gp: smpl_t = (*bt).gp;
    let mut bp: smpl_t = (*bt).bp;
    let mut rp: smpl_t = (*bt).rp;
    let mut rp1: smpl_t = (*bt).rp1;
    let mut rp2: smpl_t = (*bt).rp2;
    let mut laglen: uint_t = (*(*bt).rwv).length;
    let mut acflen: uint_t = (*(*bt).acf).length;
    let mut step: uint_t = (*bt).step;
    let mut acf: *mut fvec_t = (*bt).acf;
    let mut acfout: *mut fvec_t = (*bt).acfout;
    if gp != 0. {
        // compute shift invariant comb filterbank
        fvec_zeros(acfout);
        i = 1 as libc::c_int as uint_t;
        while i < laglen.wrapping_sub(1 as libc::c_int as libc::c_uint) {
            a = 1 as libc::c_int as uint_t;
            while a <= (*bt).timesig {
                b = 1 as libc::c_int as uint_t;
                while b < (2 as libc::c_int as libc::c_uint).wrapping_mul(a) {
                    let ref mut fresh2 = *(*acfout).data.offset(i as isize);
                    *fresh2 +=
                        *(*acf).data.offset(i.wrapping_mul(a).wrapping_add(b).wrapping_sub(1
                                                                                               as
                                                                                               libc::c_int
                                                                                               as
                                                                                               libc::c_uint)
                                                as isize);
                    b = b.wrapping_add(1)
                }
                a = a.wrapping_add(1)
            }
            i = i.wrapping_add(1)
        }
        // since gp is set, gwv has been computed in previous checkstate
        fvec_weight(acfout, (*bt).gwv);
        gp = fvec_quadratic_peak_pos(acfout, fvec_max_elem(acfout))
    } else {
        //still only using general model
        gp = 0 as libc::c_int as smpl_t
    }
    //now look for step change - i.e. a difference between gp and rp that
  // is greater than 2*constthresh - always true in first case, since gp = 0
    if counter == 0 as libc::c_int {
        if fabsf(gp - rp) as libc::c_double >
               2.0f64 * (*bt).g_var as libc::c_double {
            flagstep = 1 as libc::c_int as uint_t;
            counter = 3 as libc::c_int // have observed  step change.
            // setup 3 frame counter
        } else { flagstep = 0 as libc::c_int as uint_t }
    }
    //i.e. 3rd frame after flagstep initially set
    if counter == 1 as libc::c_int &&
           flagstep == 1 as libc::c_int as libc::c_uint {
        //check for consistency between previous beatperiod values
        if fabsf(2 as libc::c_int as libc::c_float * rp - rp1 - rp2) <
               (*bt).g_var {
            //if true, can activate context dependent model
            flagconst = 1 as libc::c_int as uint_t;
            counter = 0 as libc::c_int
            // reset counter and flagstep
        } else {
            //if not consistent, then don't flag consistency!
            flagconst = 0 as libc::c_int as uint_t;
            counter = 2 as libc::c_int
            // let it look next time
        }
    } else if counter > 0 as libc::c_int {
        //if counter doesn't = 1,
        counter = counter - 1 as libc::c_int
    }
    rp2 = rp1;
    rp1 = rp;
    if flagconst != 0 {
        /* first run of new hypothesis */
        gp = rp;
        (*bt).timesig = fvec_gettimesig(acf, acflen, gp as uint_t);
        j = 0 as libc::c_int as uint_t;
        while j < laglen {
            *(*(*bt).gwv).data.offset(j as isize) =
                expf((-0.5f64 *
                          ((j as libc::c_double + 1.0f64 -
                                gp as libc::c_double) as smpl_t *
                               (j as libc::c_double + 1.0f64 -
                                    gp as libc::c_double) as smpl_t) as
                              libc::c_double /
                          ((*bt).g_var * (*bt).g_var) as libc::c_double) as
                         libc::c_float);
            j = j.wrapping_add(1)
        }
        flagconst = 0 as libc::c_int as uint_t;
        bp = gp;
        /* flat phase weighting */
        fvec_ones((*bt).phwv);
    } else if (*bt).timesig != 0 {
        /* context dependant model */
        bp = gp;
        /* gaussian phase weighting */
        if step as libc::c_float > (*bt).lastbeat {
            j = 0 as libc::c_int as uint_t;
            while j < (2 as libc::c_int as libc::c_uint).wrapping_mul(laglen)
                  {
                *(*(*bt).phwv).data.offset(j as isize) =
                    expf((-0.5f64 *
                              ((1.0f64 + j as libc::c_double -
                                    step as libc::c_double +
                                    (*bt).lastbeat as libc::c_double) as
                                   smpl_t *
                                   (1.0f64 + j as libc::c_double -
                                        step as libc::c_double +
                                        (*bt).lastbeat as libc::c_double) as
                                       smpl_t) as libc::c_double /
                              (bp as libc::c_double / 8.0f64)) as
                             libc::c_float);
                j = j.wrapping_add(1)
            }
        } else {
            //AUBIO_DBG("NOT using phase weighting as step is %d and lastbeat %d \n",
      //                step,bt->lastbeat);
            fvec_ones((*bt).phwv);
        }
    } else {
        /* initial state */
        bp = rp;
        /* flat phase weighting */
        fvec_ones((*bt).phwv);
    }
    /* do some further checks on the final bp value */
    /* if tempo is > 206 bpm, half it */
    while (0 as libc::c_int as libc::c_float) < bp &&
              bp < 25 as libc::c_int as libc::c_float {
        /* AUBIO_BEAT_WARNINGS */
        bp = bp * 2 as libc::c_int as libc::c_float
    }
    //AUBIO_DBG("tempo:\t%3.5f bpm | ", 5168./bp);
    /* smoothing */
  //bp = (uint_t) (0.8 * (smpl_t)bp + 0.2 * (smpl_t)bp2);
  //AUBIO_DBG("tempo:\t%3.5f bpm smoothed | bp2 %d | bp %d | ", 5168./bp, bp2, bp);
  //bp2 = bp;
  //AUBIO_DBG("time signature: %d \n", bt->timesig);
    (*bt).counter = counter;
    (*bt).flagstep = flagstep;
    (*bt).gp = gp;
    (*bt).bp = bp;
    (*bt).rp1 = rp1;
    (*bt).rp2 = rp2;
}
/* * get current beat period in samples

  \param bt beat tracking object

  Returns the currently observed period, in samples, or 0 if no consistent
  value is found.

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_beattracking_get_period(mut bt:
                                                           *const aubio_beattracking_t)
 -> smpl_t {
    return (*bt).hop_size as libc::c_float * (*bt).bp;
}
/* * get current beat period in seconds

  \param bt beat tracking object

  Returns the currently observed period, in seconds, or 0 if no consistent
  value is found.

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_beattracking_get_period_s(mut bt:
                                                             *const aubio_beattracking_t)
 -> smpl_t {
    return aubio_beattracking_get_period(bt) / (*bt).samplerate as smpl_t;
}
/* * get current tempo in bpm

  \param bt beat tracking object

  Returns the currently observed tempo, in beats per minutes, or 0 if no
  consistent value is found.

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_beattracking_get_bpm(mut bt:
                                                        *const aubio_beattracking_t)
 -> smpl_t {
    if (*bt).bp != 0 as libc::c_int as libc::c_float {
        return (60.0f64 /
                    aubio_beattracking_get_period_s(bt) as libc::c_double) as
                   smpl_t
    } else { return 0.0f64 as smpl_t };
}
/* * get current tempo confidence

  \param bt beat tracking object

  Returns the confidence with which the tempo has been observed, 0 if no
  consistent value is found.

*/
#[no_mangle]
pub unsafe extern "C" fn aubio_beattracking_get_confidence(mut bt:
                                                               *const aubio_beattracking_t)
 -> smpl_t {
    if (*bt).gp != 0. {
        let mut acf_sum: smpl_t = fvec_sum((*bt).acfout);
        if acf_sum as libc::c_double != 0.0f64 {
            return fvec_quadratic_peak_mag((*bt).acfout, (*bt).gp) / acf_sum
        }
    }
    return 0.0f64 as smpl_t;
}
