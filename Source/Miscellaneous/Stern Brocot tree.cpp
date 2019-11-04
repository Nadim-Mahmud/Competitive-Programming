/**
    Name : Stern Brocot tree
    here n = numerator d = denominator
    prn = precissions numerator prd =
    O(log n)
    stern brocot tree can generate any fractional number
*/

void stern_brocot(ll n, ll d,ll prn,ll prd){
    ll ln=0,ld=1,rn=1,rd=0,mn,md;
    while(1){
        mn = ln+rn;
        md = ld+rd;
        if(n*md<d*mn){
            rn = mn;
            rd = md;
            /** condition for precissition and breaking condition
                it could be different for different problem
                here we are finding a/b which is near to c/d and a/b<c/d
                and b is as small as possible
            */
            if(prd*(d*mn-n*md)<=d*md*prn&&mn){
                x = mn;
                y = md;
                break;
            }
        }
        else{
            ln = mn;
            ld = md;
        }
    }
}
