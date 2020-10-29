from constants import PI

# To range ( -PI, PI ]
def toPI( x ) :
    while x <= -PI or PI < x :
        if ( x <= -PI ) :
            x += 2*PI
        else :
            x -= 2*PI
    return x



# To range [0, 360)
def to360( x ) :
    while x < 0 or 360 <= x :
        if x < 0 :
            x += 360
        else :
            x -= 360
    return x