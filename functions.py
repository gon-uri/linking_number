from scipy import interpolate
import numpy as np

def _rect_inter_inner(x1,x2):
    n1=x1.shape[0]-1
    n2=x2.shape[0]-1
    X1=np.c_[x1[:-1],x1[1:]]
    X2=np.c_[x2[:-1],x2[1:]]
    S1=np.tile(X1.min(axis=1),(n2,1)).T
    S2=np.tile(X2.max(axis=1),(n1,1))
    S3=np.tile(X1.max(axis=1),(n2,1)).T
    S4=np.tile(X2.min(axis=1),(n1,1))
    return S1,S2,S3,S4

def _rectangle_intersection_(x1,y1,x2,y2):
    S1,S2,S3,S4=_rect_inter_inner(x1,x2)
    S5,S6,S7,S8=_rect_inter_inner(y1,y2)

    C1=np.less_equal(S1,S2)
    C2=np.greater_equal(S3,S4)
    C3=np.less_equal(S5,S6)
    C4=np.greater_equal(S7,S8)

    ii,jj=np.nonzero(C1 & C2 & C3 & C4)
    return ii,jj

def intersection(x1,y1,x2,y2):
    """
    INTERSECTIONS Intersections of curves.
    Computes the (x,y) locations where two curves intersect.  The curves
    can be broken with NaNs or have vertical segments.

    This part of the code was taken from: https://github.com/sukhbinder/intersection
    """
    ii,jj=_rectangle_intersection_(x1,y1,x2,y2)
    n=len(ii)

    dxy1=np.diff(np.c_[x1,y1],axis=0)
    dxy2=np.diff(np.c_[x2,y2],axis=0)

    T=np.zeros((4,n))
    AA=np.zeros((4,4,n))
    AA[0:2,2,:]=-1
    AA[2:4,3,:]=-1
    AA[0::2,0,:]=dxy1[ii,:].T
    AA[1::2,1,:]=dxy2[jj,:].T

    BB=np.zeros((4,n))
    BB[0,:]=-x1[ii].ravel()
    BB[1,:]=-x2[jj].ravel()
    BB[2,:]=-y1[ii].ravel()
    BB[3,:]=-y2[jj].ravel()

    for i in range(n):
        try:
            T[:,i]=np.linalg.solve(AA[:,:,i],BB[:,i])
        except:
            T[:,i]=np.NaN


    in_range= (T[0,:] >=0) & (T[1,:] >=0) & (T[0,:] <=1) & (T[1,:] <=1)

    xy0=T[2:,in_range]
    xy0=xy0.T
    return xy0[:,0],xy0[:,1]

def _closer_point_index(curve,point):
    dist_vec = np.subtract(curve,point)
    normas = np.asarray([np.linalg.norm(vec) for vec in dist_vec])
    ind = np.argmin(normas)
    return ind

def _mad(a, axis=None):
    """
    Compute *Median Absolute Deviation* of an array along given axis.
    """

    # Median along given axis, but *keeping* the reduced axis so that
    # result can still broadcast against a.
    med = np.median(a, axis=axis, keepdims=True)
    mad = np.median(np.abs(a - med), axis=axis)  # MAD along given axis

    return mad


def linking_number(curve_1,curve_2, projection = 'AUTO',puntos_curva=5000, margin=10, verbose=False):
    """
    Computes the *Linking Number* between two three dimensional curves.

    Inputs
    curve_1: array_like
        array of shape (3, N1) with N1 the number of points of curve 1
    curve_2: array_like
        array of shape (3, N2) with N2 the number of points of curve 2
    projection: {'XY','ZX','YZ','AUTO'}, optional
        Projection plane where the intersections will be computed (default 'AUTO')
    puntos_curva: int, optional
        Number of interpolation for B-spline method (default 5000)
    margin:  int, optional
        Distance from intercetption to compute the vector (default 10)
    verbose: boolean, optional
        Print information about the computation (default False)

    Outputs
    total: float
        Linking Number (it should always be an integer, if not, check changin parameters)
    coords_1: 
        Coordinates of intesection points in the fist dimension of the projection plane (XY','ZX' or'YZ').
    coords_2: 
        Coordinates of intesection points in the second dimension of the projection plane (XY','ZX' or'YZ').

    """

    # We look for the coordinates where the curves are less 'collapsed'
    # For this we calculate the 'mad' (Median Absolute Deviation) in each dimension
    # We project on the dimensions with the highest value of mad
    if projection == 'AUTO':
        mad_x = _mad(curve_1[0]) + _mad(curve_2[0])
        mad_y = _mad(curve_1[1]) + _mad(curve_2[1])
        mad_z = _mad(curve_1[2]) + _mad(curve_2[2])
        if mad_z <= mad_x and mad_z <= mad_y:
            curve_1_projected = curve_1
            curve_2_projected = curve_2
            if verbose == True:
                print('Auto-Projection: XY')
        else:
            if mad_y <= mad_x and mad_y <= mad_x:
                curve_1_projected = [curve_1[2],curve_1[0],curve_1[1]]
                curve_2_projected = [curve_2[2],curve_2[0],curve_2[1]]
                if verbose == True:
                    print('Auto-Projection: ZX')

            else:
                if mad_x <= mad_y and mad_x <= mad_z:
                    curve_1_projected = [curve_1[1],curve_1[2],curve_1[0]]
                    curve_2_projected = [curve_2[1],curve_2[2],curve_2[0]]
                    if verbose == True:
                        print('Auto-Projection: YZ')

    if projection == 'XY':
        if verbose == True:
                    print('Projection: XY')
        curve_1_projected = curve_1
        curve_2_projected = curve_2

    if projection == 'ZX':
        if verbose == True:
                    print('Projection: ZX')
        curve_1_projected = [curve_1[2],curve_1[0],curve_1[1]]
        curve_2_projected = [curve_2[2],curve_2[0],curve_2[1]]

    if projection == 'YZ':
        if verbose == True:
                    print('Projection: YZ')
        curve_1_projected = [curve_1[1],curve_1[2],curve_1[0]]
        curve_2_projected = [curve_2[1],curve_2[2],curve_2[0]]
        
    # Puntos Curva:
    # Value that determines the resolution with which we are going to resample the curves
    # (the higher the resolution, the better the result of '_closer_point_index')
   
    # We define the vector with which we will interpolate the curves
    pasos_curva = np.divide(1,puntos_curva)
    unew = np.arange(0, 1+np.divide(pasos_curva,2), pasos_curva)

    # Interpolation of both curves
    tck, u = interpolate.splprep(curve_1_projected, s=0,k=1)
    [c1_dim1, c1_dim2,c1_dim3] = interpolate.splev(unew, tck)
    tck, u = interpolate.splprep(curve_2_projected, s=0,k=1)
    [c2_dim1, c2_dim2,c2_dim3] = interpolate.splev(unew, tck)

    # We compute the intersection between both curves
    coords_1,coords_2=intersection(c1_dim1,c1_dim2,c2_dim1,c2_dim2)

    sign = []
    arriba = []
    total = 0

    # We iterate over each intersection found
    for num in range(len(coords_1)):

        point = np.asarray([coords_1[num],coords_2[num]])

        if verbose == True:
            print('\nDim1:',point[0])
            print('Dim2:',point[1])
            
        curve_1 = np.array(list(zip(c1_dim1,c1_dim2)))
        curve_2 = np.array(list(zip(c2_dim1,c2_dim2)))

        inter_index_1 = _closer_point_index(curve_1,point)
        inter_index_2 = _closer_point_index(curve_2,point)

        
        # If the intersection is near the end of the curve
        # the sum of 'margin' can exceed the length of it
        if inter_index_1 + margin >= len(curve_1):
            final_1 = curve_1[inter_index_1+margin-len(curve_1)]
        else:
            final_1 = curve_1[inter_index_1+margin]
        if inter_index_2 + margin >= len(curve_2):
            final_2 = curve_2[inter_index_2+margin-len(curve_2)]
        else:
            final_2 = curve_2[inter_index_2+margin]

        # In the case of the beginning of the curve it resolves itself
        # (a negative index checks the vector backwards)
        inicio_1 = curve_1[inter_index_1-margin]
        inicio_2 = curve_2[inter_index_2-margin]

        # Definimos un vector que indique la direccion y sentido del flujo en el punto de intersección
        # (lo hacemos restando un punto de un instante anterior a un punto de un instante posterior)
        flecha_1 = np.subtract(final_1,inicio_1)
        flecha_2 = np.subtract(final_2,inicio_2)

        #print('Altura_1 : ',c1_dim3[inter_index_1])
        #print('Altura_2 : ',c2_dim3[inter_index_2])

        # We check which curve passes overhead at the intersection (In the coordinate that we do NOT use in the projection)
        # We compute the vector product to know if it goes clockwise or counter-clockwise
        if c1_dim3[inter_index_1] > c2_dim3[inter_index_2]:
            if verbose == True:
                print('Upper curve: Curve 1')
            arriba.append(0)
            prod_vec = np.cross(flecha_1,flecha_2)

        else:
            if verbose == True:
                print('Upper curve: Curve 2')
            arriba.append(1)
            prod_vec = np.cross(flecha_2,flecha_1)
            
        if verbose == True:
            print('Vector Product:',prod_vec)

        # We add to the total according to the result
        if prod_vec > 0 :
            sign.append(1)
            total = total + 0.5
        else:
            sign.append(-1)
            total = total - 0.5

    if verbose == True:
        print('\nLinking number: ', total)

    return total,coords_1,coords_2
