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
/* * set feedback and feedforward coefficients of a A-weighting filter

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
pub unsafe extern "C" fn aubio_filter_set_a_weighting(mut f:
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
                  b"AUBIO ERROR: aubio_filter: failed setting A-weighting with samplerate %d\n\x00"
                      as *const u8 as *const libc::c_char, samplerate);
        return AUBIO_FAIL as libc::c_int as uint_t
    }
    if f.is_null() {
        aubio_log(AUBIO_LOG_ERR as libc::c_int,
                  b"AUBIO ERROR: aubio_filter: failed setting A-weighting with filter NULL\n\x00"
                      as *const u8 as *const libc::c_char);
        return AUBIO_FAIL as libc::c_int as uint_t
    }
    order = aubio_filter_get_order(f);
    if order != 7 as libc::c_int as libc::c_uint {
        aubio_log(AUBIO_LOG_ERR as libc::c_int,
                  b"AUBIO ERROR: aubio_filter: order of A-weighting filter must be 7, not %d\n\x00"
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
                6.306209468238731519207362907764036208391189575195312500e-01f64;
            *b.offset(1 as libc::c_int as isize) =
                -1.261241893647746525886077506584115326404571533203125000e+00f64;
            *b.offset(2 as libc::c_int as isize) =
                -6.306209468238730408984338282607495784759521484375000000e-01f64;
            *b.offset(3 as libc::c_int as isize) =
                2.522483787295493051772155013168230652809143066406250000e+00f64;
            *b.offset(4 as libc::c_int as isize) =
                -6.306209468238730408984338282607495784759521484375000000e-01f64;
            *b.offset(5 as libc::c_int as isize) =
                -1.261241893647746525886077506584115326404571533203125000e+00f64;
            *b.offset(6 as libc::c_int as isize) =
                6.306209468238731519207362907764036208391189575195312500e-01f64;
            *a.offset(0 as libc::c_int as isize) =
                1.000000000000000000000000000000000000000000000000000000e+00f64;
            *a.offset(1 as libc::c_int as isize) =
                -2.128467193009123015201566886389628052711486816406250000e+00f64;
            *a.offset(2 as libc::c_int as isize) =
                2.948668980101234460278192273108288645744323730468750000e-01f64;
            *a.offset(3 as libc::c_int as isize) =
                1.824183830735050637628091863007284700870513916015625000e+00f64;
            *a.offset(4 as libc::c_int as isize) =
                -8.056628943119792385374466903158463537693023681640625000e-01f64;
            *a.offset(5 as libc::c_int as isize) =
                -3.947497982842933517133587884018197655677795410156250000e-01f64;
            *a.offset(6 as libc::c_int as isize) =
                2.098548546080332977137317129745497368276119232177734375e-01f64
        }
        11025 => {
            *b.offset(0 as libc::c_int as isize) =
                6.014684165832374640459079273568931967020034790039062500e-01f64;
            *b.offset(1 as libc::c_int as isize) =
                -1.202936833166475150136420779745094478130340576171875000e+00f64;
            *b.offset(2 as libc::c_int as isize) =
                -6.014684165832373530236054648412391543388366699218750000e-01f64;
            *b.offset(3 as libc::c_int as isize) =
                2.405873666332950300272841559490188956260681152343750000e+00f64;
            *b.offset(4 as libc::c_int as isize) =
                -6.014684165832373530236054648412391543388366699218750000e-01f64;
            *b.offset(5 as libc::c_int as isize) =
                -1.202936833166475150136420779745094478130340576171875000e+00f64;
            *b.offset(6 as libc::c_int as isize) =
                6.014684165832374640459079273568931967020034790039062500e-01f64;
            *a.offset(0 as libc::c_int as isize) =
                1.000000000000000000000000000000000000000000000000000000e+00f64;
            *a.offset(1 as libc::c_int as isize) =
                -2.463578747722854345170162559952586889266967773437500000e+00f64;
            *a.offset(2 as libc::c_int as isize) =
                1.096799662705210121060872552334330976009368896484375000e+00f64;
            *a.offset(3 as libc::c_int as isize) =
                1.381222210556041218865175324026495218276977539062500000e+00f64;
            *a.offset(4 as libc::c_int as isize) =
                -1.013875696476876031582037285261321812868118286132812500e+00f64;
            *a.offset(5 as libc::c_int as isize) =
                -1.839132734476921215982514468123554252088069915771484375e-01f64;
            *a.offset(6 as libc::c_int as isize) =
                1.833526393172056623281918064094497822225093841552734375e-01f64
        }
        16000 => {
            *b.offset(0 as libc::c_int as isize) =
                5.314898298235570806014038680586963891983032226562500000e-01f64;
            *b.offset(1 as libc::c_int as isize) =
                -1.062979659647114161202807736117392778396606445312500000e+00f64;
            *b.offset(2 as libc::c_int as isize) =
                -5.314898298235570806014038680586963891983032226562500000e-01f64;
            *b.offset(3 as libc::c_int as isize) =
                2.125959319294228322405615472234785556793212890625000000e+00f64;
            *b.offset(4 as libc::c_int as isize) =
                -5.314898298235570806014038680586963891983032226562500000e-01f64;
            *b.offset(5 as libc::c_int as isize) =
                -1.062979659647114161202807736117392778396606445312500000e+00f64;
            *b.offset(6 as libc::c_int as isize) =
                5.314898298235570806014038680586963891983032226562500000e-01f64;
            *a.offset(0 as libc::c_int as isize) =
                1.000000000000000000000000000000000000000000000000000000e+00f64;
            *a.offset(1 as libc::c_int as isize) =
                -2.867832572992162987191022693878039717674255371093750000e+00f64;
            *a.offset(2 as libc::c_int as isize) =
                2.221144410202312347024644623161293566226959228515625000e+00f64;
            *a.offset(3 as libc::c_int as isize) =
                4.552683347886614662058946123579517006874084472656250000e-01f64;
            *a.offset(4 as libc::c_int as isize) =
                -9.833868636162828025248927588108927011489868164062500000e-01f64;
            *a.offset(5 as libc::c_int as isize) =
                5.592994142413361402521587706360151059925556182861328125e-02f64;
            *a.offset(6 as libc::c_int as isize) =
                1.188781038285612462468421313133148942142724990844726562e-01f64
        }
        22050 => {
            *b.offset(0 as libc::c_int as isize) =
                4.492998504299193784916610638902056962251663208007812500e-01f64;
            *b.offset(1 as libc::c_int as isize) =
                -8.985997008598388680056245902960654348134994506835937500e-01f64;
            *b.offset(2 as libc::c_int as isize) =
                -4.492998504299192674693586013745516538619995117187500000e-01f64;
            *b.offset(3 as libc::c_int as isize) =
                1.797199401719677958055854105623438954353332519531250000e+00f64;
            *b.offset(4 as libc::c_int as isize) =
                -4.492998504299192674693586013745516538619995117187500000e-01f64;
            *b.offset(5 as libc::c_int as isize) =
                -8.985997008598388680056245902960654348134994506835937500e-01f64;
            *b.offset(6 as libc::c_int as isize) =
                4.492998504299193784916610638902056962251663208007812500e-01f64;
            *a.offset(0 as libc::c_int as isize) =
                1.000000000000000000000000000000000000000000000000000000e+00f64;
            *a.offset(1 as libc::c_int as isize) =
                -3.229078805225074955131958631682209670543670654296875000e+00f64;
            *a.offset(2 as libc::c_int as isize) =
                3.354494881236033787530459449044428765773773193359375000e+00f64;
            *a.offset(3 as libc::c_int as isize) =
                -7.317843680657351024265722116979304701089859008789062500e-01f64;
            *a.offset(4 as libc::c_int as isize) =
                -6.271627581807257545420952737913466989994049072265625000e-01f64;
            *a.offset(5 as libc::c_int as isize) =
                1.772142005020879151899748649157118052244186401367187500e-01f64;
            *a.offset(6 as libc::c_int as isize) =
                5.631716697383508385410522123493137769401073455810546875e-02f64
        }
        24000 => {
            *b.offset(0 as libc::c_int as isize) =
                4.256263892891054001488271296693710610270500183105468750e-01f64;
            *b.offset(1 as libc::c_int as isize) =
                -8.512527785782106892753517968230880796909332275390625000e-01f64;
            *b.offset(2 as libc::c_int as isize) =
                -4.256263892891054556599783609271980822086334228515625000e-01f64;
            *b.offset(3 as libc::c_int as isize) =
                1.702505557156421378550703593646176159381866455078125000e+00f64;
            *b.offset(4 as libc::c_int as isize) =
                -4.256263892891054556599783609271980822086334228515625000e-01f64;
            *b.offset(5 as libc::c_int as isize) =
                -8.512527785782106892753517968230880796909332275390625000e-01f64;
            *b.offset(6 as libc::c_int as isize) =
                4.256263892891054001488271296693710610270500183105468750e-01f64;
            *a.offset(0 as libc::c_int as isize) =
                1.000000000000000000000000000000000000000000000000000000e+00f64;
            *a.offset(1 as libc::c_int as isize) =
                -3.325996004241962733516402295208536088466644287109375000e+00f64;
            *a.offset(2 as libc::c_int as isize) =
                3.677161079286316969216841243905946612358093261718750000e+00f64;
            *a.offset(3 as libc::c_int as isize) =
                -1.106476076828482035807610373012721538543701171875000000e+00f64;
            *a.offset(4 as libc::c_int as isize) =
                -4.726706734908718843257702246773988008499145507812500000e-01f64;
            *a.offset(5 as libc::c_int as isize) =
                1.861941760230954034938122276798821985721588134765625000e-01f64;
            *a.offset(6 as libc::c_int as isize) =
                4.178771337829546850262119050967157818377017974853515625e-02f64
        }
        32000 => {
            *b.offset(0 as libc::c_int as isize) =
                3.434583386824304196416335344110848382115364074707031250e-01f64;
            *b.offset(1 as libc::c_int as isize) =
                -6.869166773648609503055695313378237187862396240234375000e-01f64;
            *b.offset(2 as libc::c_int as isize) =
                -3.434583386824303641304823031532578170299530029296875000e-01f64;
            *b.offset(3 as libc::c_int as isize) =
                1.373833354729721900611139062675647437572479248046875000e+00f64;
            *b.offset(4 as libc::c_int as isize) =
                -3.434583386824303641304823031532578170299530029296875000e-01f64;
            *b.offset(5 as libc::c_int as isize) =
                -6.869166773648609503055695313378237187862396240234375000e-01f64;
            *b.offset(6 as libc::c_int as isize) =
                3.434583386824304196416335344110848382115364074707031250e-01f64;
            *a.offset(0 as libc::c_int as isize) =
                1.000000000000000000000000000000000000000000000000000000e+00f64;
            *a.offset(1 as libc::c_int as isize) =
                -3.656446043233668063976438133977353572845458984375000000e+00f64;
            *a.offset(2 as libc::c_int as isize) =
                4.831468450652579349480220116674900054931640625000000000e+00f64;
            *a.offset(3 as libc::c_int as isize) =
                -2.557597496581567764195597192156128585338592529296875000e+00f64;
            *a.offset(4 as libc::c_int as isize) =
                2.533680394205302666144064005493419244885444641113281250e-01f64;
            *a.offset(5 as libc::c_int as isize) =
                1.224430322452567110325105659285327419638633728027343750e-01f64;
            *a.offset(6 as libc::c_int as isize) =
                6.764072168342137418572956875095769646577537059783935547e-03f64
        }
        44100 => {
            *b.offset(0 as libc::c_int as isize) =
                2.557411252042575133813784304948057979345321655273437500e-01f64;
            *b.offset(1 as libc::c_int as isize) =
                -5.114822504085150267627568609896115958690643310546875000e-01f64;
            *b.offset(2 as libc::c_int as isize) =
                -2.557411252042575133813784304948057979345321655273437500e-01f64;
            *b.offset(3 as libc::c_int as isize) =
                1.022964500817030053525513721979223191738128662109375000e+00f64;
            *b.offset(4 as libc::c_int as isize) =
                -2.557411252042575133813784304948057979345321655273437500e-01f64;
            *b.offset(5 as libc::c_int as isize) =
                -5.114822504085150267627568609896115958690643310546875000e-01f64;
            *b.offset(6 as libc::c_int as isize) =
                2.557411252042575133813784304948057979345321655273437500e-01f64;
            *a.offset(0 as libc::c_int as isize) =
                1.000000000000000000000000000000000000000000000000000000e+00f64;
            *a.offset(1 as libc::c_int as isize) =
                -4.019576181115832369528106937650591135025024414062500000e+00f64;
            *a.offset(2 as libc::c_int as isize) =
                6.189406442920693862674852425698190927505493164062500000e+00f64;
            *a.offset(3 as libc::c_int as isize) =
                -4.453198903544116404873420833609998226165771484375000000e+00f64;
            *a.offset(4 as libc::c_int as isize) =
                1.420842949621876627475103305187076330184936523437500000e+00f64;
            *a.offset(5 as libc::c_int as isize) =
                -1.418254738303044160119270600262098014354705810546875000e-01f64;
            *a.offset(6 as libc::c_int as isize) =
                4.351177233495117681327801761881346465088427066802978516e-03f64
        }
        48000 => {
            *b.offset(0 as libc::c_int as isize) =
                2.343017922995132285013397677175817079842090606689453125e-01f64;
            *b.offset(1 as libc::c_int as isize) =
                -4.686035845990265125138307666929904371500015258789062500e-01f64;
            *b.offset(2 as libc::c_int as isize) =
                -2.343017922995132007457641520886681973934173583984375000e-01f64;
            *b.offset(3 as libc::c_int as isize) =
                9.372071691980530250276615333859808743000030517578125000e-01f64;
            *b.offset(4 as libc::c_int as isize) =
                -2.343017922995132007457641520886681973934173583984375000e-01f64;
            *b.offset(5 as libc::c_int as isize) =
                -4.686035845990265125138307666929904371500015258789062500e-01f64;
            *b.offset(6 as libc::c_int as isize) =
                2.343017922995132285013397677175817079842090606689453125e-01f64;
            *a.offset(0 as libc::c_int as isize) =
                1.000000000000000000000000000000000000000000000000000000e+00f64;
            *a.offset(1 as libc::c_int as isize) =
                -4.113043408775872045168853219365701079368591308593750000e+00f64;
            *a.offset(2 as libc::c_int as isize) =
                6.553121752655050258340452273841947317123413085937500000e+00f64;
            *a.offset(3 as libc::c_int as isize) =
                -4.990849294163385074796224216697737574577331542968750000e+00f64;
            *a.offset(4 as libc::c_int as isize) =
                1.785737302937575599059982778271660208702087402343750000e+00f64;
            *a.offset(5 as libc::c_int as isize) =
                -2.461905953194876706113802811159985139966011047363281250e-01f64;
            *a.offset(6 as libc::c_int as isize) =
                1.122425003323123879339640041052916785702109336853027344e-02f64
        }
        88200 => {
            *b.offset(0 as libc::c_int as isize) =
                1.118876366882113199130444058937428053468465805053710938e-01f64;
            *b.offset(1 as libc::c_int as isize) =
                -2.237752733764226120705131961585721001029014587402343750e-01f64;
            *b.offset(2 as libc::c_int as isize) =
                -1.118876366882113337908322137081995606422424316406250000e-01f64;
            *b.offset(3 as libc::c_int as isize) =
                4.475505467528452241410263923171442002058029174804687500e-01f64;
            *b.offset(4 as libc::c_int as isize) =
                -1.118876366882113337908322137081995606422424316406250000e-01f64;
            *b.offset(5 as libc::c_int as isize) =
                -2.237752733764226120705131961585721001029014587402343750e-01f64;
            *b.offset(6 as libc::c_int as isize) =
                1.118876366882113199130444058937428053468465805053710938e-01f64;
            *a.offset(0 as libc::c_int as isize) =
                1.000000000000000000000000000000000000000000000000000000e+00f64;
            *a.offset(1 as libc::c_int as isize) =
                -4.726938565651158441482948546763509511947631835937500000e+00f64;
            *a.offset(2 as libc::c_int as isize) =
                9.076897983832765248735086061060428619384765625000000000e+00f64;
            *a.offset(3 as libc::c_int as isize) =
                -9.014855113464800950850985827855765819549560546875000000e+00f64;
            *a.offset(4 as libc::c_int as isize) =
                4.852772261031594425162438710685819387435913085937500000e+00f64;
            *a.offset(5 as libc::c_int as isize) =
                -1.333877820398965186043938047077972441911697387695312500e+00f64;
            *a.offset(6 as libc::c_int as isize) =
                1.460012549591642450064199465487035922706127166748046875e-01f64
        }
        96000 => {
            *b.offset(0 as libc::c_int as isize) =
                9.951898975972744976203898659150581806898117065429687500e-02f64;
            *b.offset(1 as libc::c_int as isize) =
                -1.990379795194548995240779731830116361379623413085937500e-01f64;
            *b.offset(2 as libc::c_int as isize) =
                -9.951898975972744976203898659150581806898117065429687500e-02f64;
            *b.offset(3 as libc::c_int as isize) =
                3.980759590389097990481559463660232722759246826171875000e-01f64;
            *b.offset(4 as libc::c_int as isize) =
                -9.951898975972744976203898659150581806898117065429687500e-02f64;
            *b.offset(5 as libc::c_int as isize) =
                -1.990379795194548995240779731830116361379623413085937500e-01f64;
            *b.offset(6 as libc::c_int as isize) =
                9.951898975972744976203898659150581806898117065429687500e-02f64;
            *a.offset(0 as libc::c_int as isize) =
                1.000000000000000000000000000000000000000000000000000000e+00f64;
            *a.offset(1 as libc::c_int as isize) =
                -4.802203044225376693532325589330866932868957519531250000e+00f64;
            *a.offset(2 as libc::c_int as isize) =
                9.401807218627226347962277941405773162841796875000000000e+00f64;
            *a.offset(3 as libc::c_int as isize) =
                -9.566143943569164420637207513209432363510131835937500000e+00f64;
            *a.offset(4 as libc::c_int as isize) =
                5.309775930392619081032989925006404519081115722656250000e+00f64;
            *a.offset(5 as libc::c_int as isize) =
                -1.517333360452622237346531619550660252571105957031250000e+00f64;
            *a.offset(6 as libc::c_int as isize) =
                1.740971994228911745583587844521389342844486236572265625e-01f64
        }
        192000 => {
            *b.offset(0 as libc::c_int as isize) =
                3.433213424548713782469278044118254911154508590698242188e-02f64;
            *b.offset(1 as libc::c_int as isize) =
                -6.866426849097426177159775306790834292769432067871093750e-02f64;
            *b.offset(2 as libc::c_int as isize) =
                -3.433213424548714476358668434841092675924301147460937500e-02f64;
            *b.offset(3 as libc::c_int as isize) =
                1.373285369819485235431955061358166858553886413574218750e-01f64;
            *b.offset(4 as libc::c_int as isize) =
                -3.433213424548714476358668434841092675924301147460937500e-02f64;
            *b.offset(5 as libc::c_int as isize) =
                -6.866426849097426177159775306790834292769432067871093750e-02f64;
            *b.offset(6 as libc::c_int as isize) =
                3.433213424548713782469278044118254911154508590698242188e-02f64;
            *a.offset(0 as libc::c_int as isize) =
                1.000000000000000000000000000000000000000000000000000000e+00f64;
            *a.offset(1 as libc::c_int as isize) =
                -5.305923689674640009172890131594613194465637207031250000e+00f64;
            *a.offset(2 as libc::c_int as isize) =
                1.165952437466175695135461864992976188659667968750000000e+01f64;
            *a.offset(3 as libc::c_int as isize) =
                -1.357560092700591525272102444432675838470458984375000000e+01f64;
            *a.offset(4 as libc::c_int as isize) =
                8.828906932824192921316353022120893001556396484375000000e+00f64;
            *a.offset(5 as libc::c_int as isize) =
                -3.039490120988216581565666274400427937507629394531250000e+00f64;
            *a.offset(6 as libc::c_int as isize) =
                4.325834301870381537469256727490574121475219726562500000e-01f64
        }
        _ => {
            aubio_log(AUBIO_LOG_ERR as libc::c_int,
                      b"AUBIO ERROR: sampling rate of A-weighting filter is %d, should be one of 8000, 11025, 16000, 22050, 24000, 32000, 44100, 48000, 88200, 96000, 192000.\n\x00"
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

  A-weighting filter coefficients

  This file creates an A-weighting digital filter, which reduces low and high
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
  1.0.9-1+b1 on i386), with <pre> [b, a] = adsign(1/Fs) </pre> for various
  sampling frequencies (8000, 11025, 16000, 22050, 24000, 32000, 44100, 48000,
  88200, 96000, and 192000 Hz).

  The sampling frequency should normally be higher than 20kHz, but most common
  file sampling rates have been included for completeness.

  \example temporal/test-a_weighting.c

*/
/* * create new A-design filter

  \param samplerate sampling frequency of the signal to filter. Should be one of
  8000, 11025, 16000, 22050, 24000, 32000, 44100, 48000, 88200, 96000, and
  192000 Hz

  \return a new filter object

*/
#[no_mangle]
pub unsafe extern "C" fn new_aubio_filter_a_weighting(mut samplerate: uint_t)
 -> *mut aubio_filter_t {
    let mut f: *mut aubio_filter_t =
        new_aubio_filter(7 as libc::c_int as uint_t);
    if aubio_filter_set_a_weighting(f, samplerate) !=
           AUBIO_OK as libc::c_int as libc::c_uint {
        del_aubio_filter(f);
        return 0 as *mut aubio_filter_t
    }
    return f;
}
