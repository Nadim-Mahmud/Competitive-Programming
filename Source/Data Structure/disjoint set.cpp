
///  *** Disjoint Set Union Find [n||1]

int parent(int n)
{
    if(rp[n]==n)return n;
    return rp[n]=parent(rp[n]);
}
void setUp(int a,int b){
    rp[parent(b)]=parent(a);
}
