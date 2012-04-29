
STATE_VARS = [ 
    'X  ', 
    'Y  ', 
    'Z  ', 
    'dX ', 
    'dY ',
    'dZ ',
    'mu ',
    'J_2',
    'C_d',
    'X_1',
    'Y_1',
    'Z_1',
    'X_2',
    'Y_2',
    'Z_2',
    'X_3',
    'Y_3',
    'Z_3' ]

MEASUREMENT_VARS = [
    'RANGE',
    'RANGE_RATE',
   ]

FIRST_CALL_TO_A = [
    [ 0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ],
    [ 0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0 ],
    [ 0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0 ],
    [ -1.04419157726501e-06, 2.48803937865918e-07, 2.32119626251703e-07, 
      -3.40563941046214e-12,-6.41345925457424e-13, 7.45142995121168e-13,
      -2.05352107782618e-15, 1.25613807792021e-00,-3.95044474364429e-09,
       0,0,0,0,0,0,0,0,0 ],
    [  2.48803890285852e-07, 6.34645016698383e-07, 1.59993315478172e-06,
      -6.41345925457424e-13,-4.18877584000629e-12, 1.32798646551291e-12,
      -1.41543269839020e-14, 8.65819654046806e-00,-7.04044349428954e-09,
       0,0,0,0,0,0,0,0,0],
    [  2.32119546345998e-07, 1.59993291293434e-06, 4.09546560812063e-07,
       7.45142995121168e-13, 1.32798646551291e-12,-4.58868668542299e-12,
      -1.31824153019035e-14,-4.42495661149975e-00, 8.17988692853170e-09,
       0,0,0,0,0,0,0,0,0],
    [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ],
    [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ],
    [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ],
    [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ],
    [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ],
    [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ],
    [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ],
    [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ],
    [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ],
    [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ],
    [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ],
    [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ],
    ]

FIRST_CALL_TO_HTILDA = [
    [ -0.815628974909, 0.521493329427, 0.250587475051, 
       0, 0, 0, 
       0, 0, 0, 
       0, 0, 0, 
       0.815628974909,-0.521493329427,-0.250587475051, 
       0, 0, 0 ],
    [  0.00041849970161, 0.00129966411126, -0.00134254787406, 
      -0.815628974909, 0.521493329427, 0.250587475051, 
       0, 0, 0, 
       0, 0, 0, 
      -0.000456527599369, -0.00135914072106, 0.00134254787406,
      0, 0, 0 ],
    ]

# Position and velocity integrated for first 120 seconds
P1_POS_VEL_T_20_120 = {
     0  :  [ 757700.000000, 5222607.000000, 4851500.000000, 2213.210000, 4678.340000, -5371.300000 ],
     20 :  [ 801797.307176, 5315038.691050, 4743030.867004, 2196.361728, 4564.494933, -5475.222907 ],
     40 :  [ 845548.133491, 5405173.594176, 4632506.851200, 2178.563145, 4448.669469, -5576.780775 ],
     60 :  [ 888933.547992, 5492972.600125, 4519975.700317, 2159.821848, 4330.913624, -5675.929080 ],
     80 :  [ 931934.775863, 5578397.608745, 4405486.042278, 2140.145853, 4211.278289, -5772.624330 ],
     100 : [ 974533.206751, 5661411.546278, 4289087.364309, 2119.543589, 4089.815209, -5866.824090 ],
     120 : [1016710.403030, 5741978.382218, 4170829.991621, 2098.023900, 3966.576961, -5958.487000 ],
    }

# The values of PHI at t=20 and t=18340
P1_PHI_T_20_18340 = {
   20 : [
     [ 0.999791449912,5.10374043053e-005,4.69769058557e-005,19.9986105457,0.000344485775519,0.000315004905618,-4.18727914e-013,249.325398126,-7.88656246882e-007, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
     [ 5.10372857847e-005,1.00013105813,0.000319525139904,0.000344485711883,20.000887281,0.00212851097752,-2.84794296036e-012,1696.0678591,-1.39548873909e-006 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
     [ 4.69766877245e-005,0.000319524398352,1.00007753871,0.000315004799188,0.00212851065534,20.0005021709,-2.61717607214e-012,-920.002678521,1.64567424183e-006 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
     [ -2.08400149199e-5,5.16791577017e-6,4.72475869704e-6,0.999791733608,5.23170599968e-005,4.75146775342e-005,-4.22748335752e-014,24.8276433094,-7.87908466515e-008 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
     [ 5.16789302185e-006,1.33129491365e-005,3.19278168585e-005,5.23169414888e-005,1.00013516,0.000319007538445,-2.85653500556e-013,167.806773084,-1.38913736664e-007 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
     [ 4.72471667384e-006,3.1927673402e-005,7.53641447024e-006,4.7514459469e-005,0.000319006797104,1.00007315288,-2.60756681803e-013,-93.7197084849,1.65045181034e-007 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
     [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
     [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
     [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
     [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
     [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
     [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
     [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
     [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
     [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
     [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
     [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
     [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
     ],
    

    18340 : [
     [ -0.616427764975,-11.0902065782,-10.322380477,-4174.38626236,-9188.75335432,10552.2303097,1.82823858434e-007,-55451122.9669,1.36087787207, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
     [ -2.85151970833,-18.6536459142,-18.3027844025,-7704.60804153,-16118.7137877,18729.9135645,3.24100644125e-007,-58773945.366,2.21201761337, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
     [ 5.01767677069,34.4992296182,33.1013477828,13508.2045768,28583.4879404,-32615.1552117,-5.68520848896e-007,109778994.208,-4.57747754265, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
     [ 0.000799903195927,0.00680753691194,0.00633657817886,3.61941394161,5.59882083469,-6.4225819708,-1.11736420942e-010,9732.02290012,-0.000820498938093, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
     [ 0.00525838818842,0.0359487426774,0.0336127215768,14.0365813459,30.6797928483,-34.0340374188,-5.9283119321e-010,121949.210357,-0.00448453330932, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
     [ 0.00333619090095,0.0228874271035,0.0210626833205,8.86261097902,18.7709787521,-20.5051552194,-3.74818578849e-010,70925.1102015,-0.00303132513226, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
     [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
     [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
     [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
     [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
     [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
     [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
     [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
     [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
     [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
     [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
     [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
     [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
    ]
    }

# Values for x_hat for three consequtive iterations
P1_XHAT = [
    [ -0.0355074965116273 ],
    [ -0.27341590598553 ],
    [ -0.179970761622627 ],
    [ 0.0409377124301569 ],
    [ 0.0327563065901537 ],
    [ -0.0147606324122311 ],
    [ -9310103.07490921 ],
    [ -6.57354368573488e-07 ],
    [ 0.108504936878731 ],
    [ 1.86320404339314e-06 ],
    [ 1.37869926339216e-06 ],
    [ -2.54491428042453e-07 ],
    [ -10.5636844179305 ],
    [ 9.98328899684054 ],
    [ 5.79442720656618 ],
    [ -5.78200041893004 ],
    [ 2.3443163443211 ],
    [ 1.5100298917132 ],
    ]

P2_XHAT = [
    [ 0.327002735286328 ],
    [ -0.148219971423151 ],
    [ -0.0809954174956853 ],
    [ -0.000317394513040911 ],
    [ -3.91734946810604e-05 ],
    [ 0.000338256138737275 ],
    [ -33324610.6637554 ],
    [ 2.99030136449693e-08 ],
    [ 0.0411509188471517 ],
    [ -7.51773924336147e-09 ],
    [ -5.56283769036147e-09 ],
    [ 1.89907559370749e-10 ],
    [ 0.555565210943936 ],
    [ 0.0200994336432552 ],
    [ 0.182189541170642 ],
    [ 0.773745083515455 ],
    [ -0.32334236251464 ],
    [ 1.46455907608498 ],
    ]

P3_XHAT = [
    [ 4.01598692561632e-08 ],
    [ -3.80136581319667e-07 ],
    [ 2.34248356314115e-07 ],
    [ -3.20528232624071e-10 ],
    [ 2.99837858339875e-10 ],
    [ 3.41973415841183e-10 ],
    [ -18.1455098503496 ],
    [ -5.76289957892531e-16 ],
    [ -1.98584649459958e-07 ],
    [ -7.51686732608076e-09 ],
    [ -5.56219243530013e-09 ],
    [ 1.90110260662791e-10 ],
    [ 5.89962559487648e-08 ],
    [ -1.60768256507461e-07 ],
    [ 2.12235549097213e-07 ],
    [ -4.67520773898634e-08 ],
    [ 3.13874466138048e-08 ],
    [ 5.60596784951685e-07 ],
    ]
