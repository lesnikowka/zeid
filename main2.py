import numpy as np

a = [#v11  v21  v31  v12  v22 v32  v13 v23 v33
    [-832., 16., 0., 400., 0., 0., 0., 0., 0.],
    [16., -832., 16., 0., 400., 0., 0., 0., 0.],
    [0., 16., -832., 0, 0., 400., 0., 0., 0.],
    [400., 0., 0., -832., 16., 0., 400., 0., 0.],
    [0., 400., 0., 16., -832., 16., 0., 400., 0.],
    [0., 0., 400., 0., 16., -832., 0., 0., 400.],
    [0., 0., 0., 400., 0., 0., -832., 16., 0.],
    [0., 0., 0., 0., 400., 0., 16., -832., 16.],
    [0., 0., 0., 0., 0., 400., 0., 16., -832.]
]

sup_val = np.linalg.eig(a)

val = [abs(val_) for val_ in sup_val[0]]

print(min(val))

a = np.array(a)

a = -a

det = np.linalg.det(a)


print(det)