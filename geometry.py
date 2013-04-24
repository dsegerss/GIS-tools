import numpy as np

SMALL_NUM=0.00000001


class Segment:
    """
    Segment defined by two points, p0 and p1 
    """
    def __init__(self,p0,p1):
        """
        @param p1: starting point as numpy array
        @param p2: end point as numpy array
        """
        self.P0=p0
        self.P1=p1

    def to2D(self):
        return Segment(self.P0[:2],self.P1[:2])

def perp2D(u,v):
    """
    Perpendicular product of two vectors
    (u.x * v.y - u.y * v.x)
    """
    return u[0] * v[1] - u[1] * v[0]


def intersect2D_2Segments(seg1,seg2):
    """
    Follows: http://geomalgorithms.com/a05-_intersect-1.html by Dan Sunday
    Check intersection of two segments
    Returns (I0,I1) where I0 is intersection point or None,
    and I1 is end point of intersection or None
    @param S1: segment 1
    @param S2: segment 2
    """

    S1=seg1.to2D()
    S2=seg2.to2D()
    
    u = S1.P1 - S1.P0
    v = S2.P1 - S2.P0
    w = S1.P0 - S2.P0
    D = perp2D(u,v)

    intersection=[None,None]
    # test if  they are parallel (includes either being a point)
    if abs(D) < SMALL_NUM:           # S1 and S2 are parallel
        if perp2D(u,w) != 0 or perp2D(v,w) != 0:
            return intersection       # they are NOT collinear

        # they are collinear or degenerate
        # check if they are degenerate  points
        du = np.dot(u,u);
        dv = np.dot(v,v);
        if du==0 and dv==0:             # both segments are points
            if S1.P0 ==  S2.P0:         
                intersection[0] = S1.P0 # they are the same point
                
        elif du==0:                     # S1 is a single point
            if inSegment(S1.P0, S2):  
                intersection[0]=S1.P0             
            
        elif (dv==0):                     # S2 a single point
            if  inSegment(S2.P0, S1): 
                intersection[0] = S2.P0
        else:
            # they are collinear segments - get  overlap (or not)
            #t0, t1 are endpoints of S1 in eqn for S2
            w2 = S1.P1 - S2.P0
            if v[0] != 0:
                t0 = w[0]/v[0]
                t1 = w2[0]/v[0]
            else:
                t0 = w[1]/ v[1]
                t1 = w2[1]/v[1]

            if t0 > t1:                  # must have t0 smaller than t1
                t=t0; t0=t1; t1=t;       # swap if not

            if t0 <= 1 and t1 >= 0:                
                t0=max(0,t0)                 # clip to min 0
                t1=min(1,t1)                 # clip to max 1
                if t0 == t1:               # intersect is a point
                    intersection[0] = S2.P0 +  t0 * v
            else:                
                # they overlap in a valid subsegment
                intersection[0] = S2.P0 + t0 * v
                intersection[1] = S2.P0 + t1 * v
    else:
        # the segments are skew and may intersect in a point
        # get the intersect parameter for S1
        sI = perp2D(v,w) / D
        if sI < 0 or sI > 1:                # no intersect with S1
            return intersection

        # get the intersect parameter for S2
        tI = perp2D(u,w) / D
        if tI < 0 or tI > 1:                # no intersect with S2
            return intersection

        intersection[0] = S1.P0 + sI * u;   # compute S1 intersect point
        return intersection



def inSegment( P, S):
    """
    Determine if a point is inside a segment
    @param P: A point
    @param S: A collinear segment S
    """
    if S.P0[0] != S.P1[0]:    # S is not  vertical
        if S.P0[0] <= P[0] and P[0] <= S.P1[0]:
            return True
        if S.P0[0] >= P[0] and P[0] >= S.P1[0]:
            return True
    else:    # S is vertical, so test y  coordinate
        if S.P0[1] <= P[1] and P[1] <= S.P1[1]:
            return True
        if S.P0[1] >= P[1] and P[1] >= S.P1[1]:
            return True
    return False