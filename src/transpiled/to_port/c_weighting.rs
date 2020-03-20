use ::libc;
extern "C" {
    pub type _aubio_filter_t;
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
    /* * @file
 * Private include file
 *
 * This file is for inclusion from _within_ the library only.
 */
    /* ********************
 *
 * External includes
 *
 */
    /* must be included before fftw3.h */
    // for CHAR_BIT, in C99 standard
    // --enable-blas=true
    /* HAVE_ACCELERATE */
    /* HAVE_BLAS */
    /* ***
 *
 * SYSTEM INTERFACE
 *
 */
    /* Memory management */
    /* file interface */
    /* strings */
    /* Error reporting */
    /* Logging */
    /* * internal logging function, defined in utils/log.c */
    #[no_mangle]
    fn aubio_log(level: sint_t, fmt: *const char_t, _: ...) -> uint_t;
    /* * returns a pointer to feedback coefficients \f$ a_i \f$

  \param f filter object to get parameters from

  \return a pointer to the \f$ a_0 ... a_i ... a_P \f$ coefficients

*/
    #[no_mangle]
    fn aubio_filter_get_feedback(f: *const aubio_filter_t) -> *mut lvec_t;
    /* * returns a pointer to feedforward coefficients \f$ b_i \f$

  \param f filter object to get coefficients from

  \return a pointer to the \f$ b_0 ... b_i ... b_P \f$ coefficients

*/
    #[no_mangle]
    fn aubio_filter_get_feedforward(f: *const aubio_filter_t) -> *mut lvec_t;
    /* * get order of the filter

  \param f filter to get order from

  \return the order of the filter

*/
    #[no_mangle]
    fn aubio_filter_get_order(f: *const aubio_filter_t) -> uint_t;
    /* * get sampling rate of the filter

  \param f filter to get sampling rate from
  \param samplerate sample rate to set the filter to

  \return the sampling rate of the filter, in Hz

*/
    #[no_mangle]
    fn aubio_filter_set_samplerate(f: *mut aubio_filter_t, samplerate: uint_t)
     -> uint_t;
    /* * create new filter object

  This function creates a new ::aubio_filter_t object, given the order of the
  filter.

  \param order order of the filter (number of coefficients)

  \return the newly created filter object

*/
    #[no_mangle]
    fn new_aubio_filter(order: uint_t) -> *mut aubio_filter_t;
    /* * delete a filter object

  \param f filter object to delete

*/
    #[no_mangle]
    fn del_aubio_filter(f: *mut aubio_filter_t);
}
/* * print format for sample in single precision */
/* * long sample format (64 bits or more) */
pub type lsmp_t = libc::c_double;
/* * print format for sample in double precision */
/* * unsigned integer */
pub type uint_t = libc::c_uint;
/* * signed integer */
pub type sint_t = libc::c_int;
/* * character */
pub type char_t = libc::c_char;
pub type C2RustUnnamed = libc::c_uint;
pub const AUBIO_FAIL: C2RustUnnamed = 1;
pub const AUBIO_OK: C2RustUnnamed = 0;
pub type aubio_log_level = libc::c_uint;
pub const AUBIO_LOG_LAST_LEVEL: aubio_log_level = 5;
pub const AUBIO_LOG_WRN: aubio_log_level = 4;
pub const AUBIO_LOG_DBG: aubio_log_level = 3;
pub const AUBIO_LOG_MSG: aubio_log_level = 2;
pub const AUBIO_LOG_INF: aubio_log_level = 1;
pub const AUBIO_LOG_ERR: aubio_log_level = 0;
#[derive(Copy, Clone)]
#[repr(C)]
pub struct lvec_t {
    pub length: uint_t,
    pub data: *mut lsmp_t,
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
/* * Digital filter

*/
pub type aubio_filter_t = _aubio_filter_t;
/* * set feedback and feedforward coefficients of a C-weighting filter

  \param f filter object to get coefficients from
  \param samplerate sampling frequency of the signal to filter. Should be one of
  8000, 11025, 16000, 22050, 24000, 32000, 44100, 48000, 88200, 96000, and
  192000 Hz

*/
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
#[no_mangle]
pub unsafe extern "C" fn aubio_filter_set_c_weighting(mut f:
                                                          *mut aubio_filter_t,
                                                      mut samplerate: uint_t)
 -> uint_t {
    let mut order: uint_t = 0;
    let mut a: *mut lsmp_t = 0 as *mut lsmp_t;
    let mut b: *mut lsmp_t = 0 as *mut lsmp_t;
    let mut as_0: *mut lvec_t = 0 as *mut lvec_t;
    let mut bs: *mut lvec_t = 0 as *mut lvec_t;
    if samplerate as sint_t <= 0 as libc::c_int {
        aubio_log(AUBIO_LOG_ERR as libc::c_int,
                  b"AUBIO ERROR: aubio_filter: failed setting C-weighting with samplerate %d\n\x00"
                      as *const u8 as *const libc::c_char, samplerate);
        return AUBIO_FAIL as libc::c_int as uint_t
    }
    if f.is_null() {
        aubio_log(AUBIO_LOG_ERR as libc::c_int,
                  b"AUBIO ERROR: aubio_filter: failed setting C-weighting with filter NULL\n\x00"
                      as *const u8 as *const libc::c_char);
        return AUBIO_FAIL as libc::c_int as uint_t
    }
    order = aubio_filter_get_order(f);
    if order != 5 as libc::c_int as libc::c_uint {
        aubio_log(AUBIO_LOG_ERR as libc::c_int,
                  b"AUBIO ERROR: aubio_filter: order of C-weighting filter must be 5, not %d\n\x00"
                      as *const u8 as *const libc::c_char, order);
        return 1 as libc::c_int as uint_t
    }
    aubio_filter_set_samplerate(f, samplerate);
    bs = aubio_filter_get_feedforward(f);
    as_0 = aubio_filter_get_feedback(f);
    b = (*bs).data;
    a = (*as_0).data;
    /* select coefficients according to sampling frequency */
    match samplerate {
        8000 => {
            *b.offset(0 as libc::c_int as isize) =
                6.782173932405135552414776611840352416038513183593750000e-01f64;
            *b.offset(1 as libc::c_int as isize) =
                0.000000000000000000000000000000000000000000000000000000e+00f64;
            *b.offset(2 as libc::c_int as isize) =
                -1.356434786481027110482955322368070483207702636718750000e+00f64;
            *b.offset(3 as libc::c_int as isize) =
                0.000000000000000000000000000000000000000000000000000000e+00f64;
            *b.offset(4 as libc::c_int as isize) =
                6.782173932405135552414776611840352416038513183593750000e-01f64;
            *a.offset(0 as libc::c_int as isize) =
                1.000000000000000000000000000000000000000000000000000000e+00f64;
            *a.offset(1 as libc::c_int as isize) =
                -6.589092811505605773447769024642184376716613769531250000e-01f64;
            *a.offset(2 as libc::c_int as isize) =
                -1.179445664897062595599663836765103042125701904296875000e+00f64;
            *a.offset(3 as libc::c_int as isize) =
                4.243329729632123736848825501510873436927795410156250000e-01f64;
            *a.offset(4 as libc::c_int as isize) =
                4.147270002091348328754349950031610205769538879394531250e-01f64
        }
        11025 => {
            *b.offset(0 as libc::c_int as isize) =
                6.002357155402652244546857218665536493062973022460937500e-01f64;
            *b.offset(1 as libc::c_int as isize) =
                0.000000000000000000000000000000000000000000000000000000e+00f64;
            *b.offset(2 as libc::c_int as isize) =
                -1.200471431080530448909371443733107298612594604492187500e+00f64;
            *b.offset(3 as libc::c_int as isize) =
                0.000000000000000000000000000000000000000000000000000000e+00f64;
            *b.offset(4 as libc::c_int as isize) =
                6.002357155402652244546857218665536493062973022460937500e-01f64;
            *a.offset(0 as libc::c_int as isize) =
                1.000000000000000000000000000000000000000000000000000000e+00f64;
            *a.offset(1 as libc::c_int as isize) =
                -8.705602141280316397242700077185872942209243774414062500e-01f64;
            *a.offset(2 as libc::c_int as isize) =
                -9.037199507150940336330791069485712796449661254882812500e-01f64;
            *a.offset(3 as libc::c_int as isize) =
                4.758433040929530011275971901341108605265617370605468750e-01f64;
            *a.offset(4 as libc::c_int as isize) =
                2.987653956523212417373258631414500996470451354980468750e-01f64
        }
        16000 => {
            *b.offset(0 as libc::c_int as isize) =
                4.971057193673903418229542694461997598409652709960937500e-01f64;
            *b.offset(1 as libc::c_int as isize) =
                0.000000000000000000000000000000000000000000000000000000e+00f64;
            *b.offset(2 as libc::c_int as isize) =
                -9.942114387347806836459085388923995196819305419921875000e-01f64;
            *b.offset(3 as libc::c_int as isize) =
                0.000000000000000000000000000000000000000000000000000000e+00f64;
            *b.offset(4 as libc::c_int as isize) =
                4.971057193673903418229542694461997598409652709960937500e-01f64;
            *a.offset(0 as libc::c_int as isize) =
                1.000000000000000000000000000000000000000000000000000000e+00f64;
            *a.offset(1 as libc::c_int as isize) =
                -1.162322939286873690889478893950581550598144531250000000e+00f64;
            *a.offset(2 as libc::c_int as isize) =
                -4.771961355734982701548574368644040077924728393554687500e-01f64;
            *a.offset(3 as libc::c_int as isize) =
                4.736145114694704227886745684372726827859878540039062500e-01f64;
            *a.offset(4 as libc::c_int as isize) =
                1.660337524309875301131711466950946487486362457275390625e-01f64
        }
        22050 => {
            *b.offset(0 as libc::c_int as isize) =
                4.033381299002754549754001800465630367398262023925781250e-01f64;
            *b.offset(1 as libc::c_int as isize) =
                0.000000000000000000000000000000000000000000000000000000e+00f64;
            *b.offset(2 as libc::c_int as isize) =
                -8.066762598005509099508003600931260734796524047851562500e-01f64;
            *b.offset(3 as libc::c_int as isize) =
                0.000000000000000000000000000000000000000000000000000000e+00f64;
            *b.offset(4 as libc::c_int as isize) =
                4.033381299002754549754001800465630367398262023925781250e-01f64;
            *a.offset(0 as libc::c_int as isize) =
                1.000000000000000000000000000000000000000000000000000000e+00f64;
            *a.offset(1 as libc::c_int as isize) =
                -1.449545371157945350404361306573264300823211669921875000e+00f64;
            *a.offset(2 as libc::c_int as isize) =
                -1.030104190885922088583015465701464563608169555664062500e-02f64;
            *a.offset(3 as libc::c_int as isize) =
                3.881857047554073680828423675848171114921569824218750000e-01f64;
            *a.offset(4 as libc::c_int as isize) =
                7.171589940116777917022972133054281584918498992919921875e-02f64
        }
        24000 => {
            *b.offset(0 as libc::c_int as isize) =
                3.786678621924967069745093795063439756631851196289062500e-01f64;
            *b.offset(1 as libc::c_int as isize) =
                0.000000000000000000000000000000000000000000000000000000e+00f64;
            *b.offset(2 as libc::c_int as isize) =
                -7.573357243849934139490187590126879513263702392578125000e-01f64;
            *b.offset(3 as libc::c_int as isize) =
                0.000000000000000000000000000000000000000000000000000000e+00f64;
            *b.offset(4 as libc::c_int as isize) =
                3.786678621924967069745093795063439756631851196289062500e-01f64;
            *a.offset(0 as libc::c_int as isize) =
                1.000000000000000000000000000000000000000000000000000000e+00f64;
            *a.offset(1 as libc::c_int as isize) =
                -1.529945307555420797029910318087786436080932617187500000e+00f64;
            *a.offset(2 as libc::c_int as isize) =
                1.283553182116208835061854642844991758465766906738281250e-01f64;
            *a.offset(3 as libc::c_int as isize) =
                3.494608072385725350272878131363540887832641601562500000e-01f64;
            *a.offset(4 as libc::c_int as isize) =
                5.217291949300089520802359288609295617789030075073242188e-02f64
        }
        32000 => {
            *b.offset(0 as libc::c_int as isize) =
                2.977986488230693340462096330156782642006874084472656250e-01f64;
            *b.offset(1 as libc::c_int as isize) =
                0.000000000000000000000000000000000000000000000000000000e+00f64;
            *b.offset(2 as libc::c_int as isize) =
                -5.955972976461386680924192660313565284013748168945312500e-01f64;
            *b.offset(3 as libc::c_int as isize) =
                0.000000000000000000000000000000000000000000000000000000e+00f64;
            *b.offset(4 as libc::c_int as isize) =
                2.977986488230693340462096330156782642006874084472656250e-01f64;
            *a.offset(0 as libc::c_int as isize) =
                1.000000000000000000000000000000000000000000000000000000e+00f64;
            *a.offset(1 as libc::c_int as isize) =
                -1.812455387128179218336754274787381291389465332031250000e+00f64;
            *a.offset(2 as libc::c_int as isize) =
                6.425013281155662614452239722595550119876861572265625000e-01f64;
            *a.offset(3 as libc::c_int as isize) =
                1.619857574578579817448087396769551560282707214355468750e-01f64;
            *a.offset(4 as libc::c_int as isize) =
                7.987649713547682189807019881300220731645822525024414062e-03f64
        }
        44100 => {
            *b.offset(0 as libc::c_int as isize) =
                2.170085619492190254220531642204150557518005371093750000e-01f64;
            *b.offset(1 as libc::c_int as isize) =
                0.000000000000000000000000000000000000000000000000000000e+00f64;
            *b.offset(2 as libc::c_int as isize) =
                -4.340171238984380508441063284408301115036010742187500000e-01f64;
            *b.offset(3 as libc::c_int as isize) =
                0.000000000000000000000000000000000000000000000000000000e+00f64;
            *b.offset(4 as libc::c_int as isize) =
                2.170085619492190254220531642204150557518005371093750000e-01f64;
            *a.offset(0 as libc::c_int as isize) =
                1.000000000000000000000000000000000000000000000000000000e+00f64;
            *a.offset(1 as libc::c_int as isize) =
                -2.134674963687040794013682898366823792457580566406250000e+00f64;
            *a.offset(2 as libc::c_int as isize) =
                1.279333533236062692139967111870646476745605468750000000e+00f64;
            *a.offset(3 as libc::c_int as isize) =
                -1.495598460893957093453821016737492755055427551269531250e-01f64;
            *a.offset(4 as libc::c_int as isize) =
                4.908700174624683852664386307651511742733418941497802734e-03f64
        }
        48000 => {
            *b.offset(0 as libc::c_int as isize) =
                1.978871200263932761398422144338837824761867523193359375e-01f64;
            *b.offset(1 as libc::c_int as isize) =
                0.000000000000000000000000000000000000000000000000000000e+00f64;
            *b.offset(2 as libc::c_int as isize) =
                -3.957742400527865522796844288677675649523735046386718750e-01f64;
            *b.offset(3 as libc::c_int as isize) =
                0.000000000000000000000000000000000000000000000000000000e+00f64;
            *b.offset(4 as libc::c_int as isize) =
                1.978871200263932761398422144338837824761867523193359375e-01f64;
            *a.offset(0 as libc::c_int as isize) =
                1.000000000000000000000000000000000000000000000000000000e+00f64;
            *a.offset(1 as libc::c_int as isize) =
                -2.219172914052801814932536217384040355682373046875000000e+00f64;
            *a.offset(2 as libc::c_int as isize) =
                1.455135878947171557129536267893854528665542602539062500e+00f64;
            *a.offset(3 as libc::c_int as isize) =
                -2.484960738877830532800317087094299495220184326171875000e-01f64;
            *a.offset(4 as libc::c_int as isize) =
                1.253882314727246607977129144728678511455655097961425781e-02f64
        }
        88200 => {
            *b.offset(0 as libc::c_int as isize) =
                9.221909851156021020734954163344809785485267639160156250e-02f64;
            *b.offset(1 as libc::c_int as isize) =
                0.000000000000000000000000000000000000000000000000000000e+00f64;
            *b.offset(2 as libc::c_int as isize) =
                -1.844381970231204204146990832668961957097053527832031250e-01f64;
            *b.offset(3 as libc::c_int as isize) =
                0.000000000000000000000000000000000000000000000000000000e+00f64;
            *b.offset(4 as libc::c_int as isize) =
                9.221909851156021020734954163344809785485267639160156250e-02f64;
            *a.offset(0 as libc::c_int as isize) =
                1.000000000000000000000000000000000000000000000000000000e+00f64;
            *a.offset(1 as libc::c_int as isize) =
                -2.785795902923448696952846148633398115634918212890625000e+00f64;
            *a.offset(2 as libc::c_int as isize) =
                2.727736758747444145711824603495188057422637939453125000e+00f64;
            *a.offset(3 as libc::c_int as isize) =
                -1.097007502819661528548067508381791412830352783203125000e+00f64;
            *a.offset(4 as libc::c_int as isize) =
                1.550674356752141103132913713125162757933139801025390625e-01f64
        }
        96000 => {
            *b.offset(0 as libc::c_int as isize) =
                8.182864044979756834585771230194950476288795471191406250e-02f64;
            *b.offset(1 as libc::c_int as isize) =
                0.000000000000000000000000000000000000000000000000000000e+00f64;
            *b.offset(2 as libc::c_int as isize) =
                -1.636572808995951366917154246038990095257759094238281250e-01f64;
            *b.offset(3 as libc::c_int as isize) =
                0.000000000000000000000000000000000000000000000000000000e+00f64;
            *b.offset(4 as libc::c_int as isize) =
                8.182864044979756834585771230194950476288795471191406250e-02f64;
            *a.offset(0 as libc::c_int as isize) =
                1.000000000000000000000000000000000000000000000000000000e+00f64;
            *a.offset(1 as libc::c_int as isize) =
                -2.856378516857566829401093855267390608787536621093750000e+00f64;
            *a.offset(2 as libc::c_int as isize) =
                2.897640237559524045707348705036565661430358886718750000e+00f64;
            *a.offset(3 as libc::c_int as isize) =
                -1.225265958339703198376469117647502571344375610351562500e+00f64;
            *a.offset(4 as libc::c_int as isize) =
                1.840048283551226071530493300087982788681983947753906250e-01f64
        }
        192000 => {
            *b.offset(0 as libc::c_int as isize) =
                2.784755468532278815940728122768632601946592330932617188e-02f64;
            *b.offset(1 as libc::c_int as isize) =
                0.000000000000000000000000000000000000000000000000000000e+00f64;
            *b.offset(2 as libc::c_int as isize) =
                -5.569510937064557631881456245537265203893184661865234375e-02f64;
            *b.offset(3 as libc::c_int as isize) =
                0.000000000000000000000000000000000000000000000000000000e+00f64;
            *b.offset(4 as libc::c_int as isize) =
                2.784755468532278815940728122768632601946592330932617188e-02f64;
            *a.offset(0 as libc::c_int as isize) =
                1.000000000000000000000000000000000000000000000000000000e+00f64;
            *a.offset(1 as libc::c_int as isize) =
                -3.333298856144166322224009491037577390670776367187500000e+00f64;
            *a.offset(2 as libc::c_int as isize) =
                4.111467536240339448738723149290308356285095214843750000e+00f64;
            *a.offset(3 as libc::c_int as isize) =
                -2.222889041651291641699117462849244475364685058593750000e+00f64;
            *a.offset(4 as libc::c_int as isize) =
                4.447204118126878991112960193277103826403617858886718750e-01f64
        }
        _ => {
            aubio_log(AUBIO_LOG_ERR as libc::c_int,
                      b"AUBIO ERROR: sampling rate of C-weighting filter is %d, should be one of 8000, 11025, 16000, 22050, 24000, 32000, 44100, 48000, 88200, 96000, 192000.\n\x00"
                          as *const u8 as *const libc::c_char, samplerate);
            return 1 as libc::c_int as uint_t
        }
    }
    return 0 as libc::c_int as uint_t;
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

  C-weighting filter coefficients

  This file creates a C-weighting digital filter, which reduces low and high
  frequencies and enhance the middle ones to reflect the ability of the human
  hearing.

  The implementation is based on the following standard:

    - IEC/CD 1672: Electroacoustics-Sound Level Meters, IEC, Geneva, Nov.  1996,
  for A- and C-weighting filters.

  See also:

    - <a href="http://en.wikipedia.org/wiki/A-weighting">A-Weighting on
  Wikipedia</a>
    - <a href="http://en.wikipedia.org/wiki/Weighting_filter">Weighting filter on
  Wikipedia</a>
    - <a href="http://www.mathworks.com/matlabcentral/fileexchange/69">Christophe
  Couvreur's 'octave' toolbox</a>

  The coefficients in this file have been computed using Christophe Couvreur's
  scripts in octave 3.0 (debian package 1:3.0.5-6+b2 with octave-signal
  1.0.9-1+b1 on i386), with <pre> [b, a] = cdsign(1/Fs) </pre> for various
  sampling frequencies (8000, 11025, 16000, 22050, 24000, 32000, 44100, 48000,
  88200, 96000, and 192000 Hz).

  The sampling frequency should normally be higher than 20kHz, but most common
  file sampling rates have been included for completeness.

  \example temporal/test-c_weighting.c

*/
/* * create new C-design filter

  \param samplerate sampling frequency of the signal to filter. Should be one of
  8000, 11025, 16000, 22050, 24000, 32000, 44100, 48000, 88200, 96000, and
  192000 Hz

  \return a new filter object

*/
#[no_mangle]
pub unsafe extern "C" fn new_aubio_filter_c_weighting(mut samplerate: uint_t)
 -> *mut aubio_filter_t {
    let mut f: *mut aubio_filter_t =
        new_aubio_filter(5 as libc::c_int as uint_t);
    if aubio_filter_set_c_weighting(f, samplerate) !=
           AUBIO_OK as libc::c_int as libc::c_uint {
        del_aubio_filter(f);
        return 0 as *mut aubio_filter_t
    }
    return f;
}
