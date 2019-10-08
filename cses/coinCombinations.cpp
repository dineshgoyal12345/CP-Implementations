#include<bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;
#define loop(i, a, b) for(lli i = a;i <= b; ++i)
#define bbkj ios::sync_with_stdio(false);cin.tie(0);cout.tie(0)
#define lli long long int 
#define p (lli)(1e9+7)
#define ordered_set tree<int, null_type, less_equal<int>, rb_tree_tag, tree_order_statistics_node_update>
using namespace std;
 
signed main(){
 
    // #ifndef ONLINE_JUDGE
    // // for getting input from input.txt
    // freopen("input.txt", "r", stdin);
    // // for writing output to output.txt
    // freopen("output.txt", "w", stdout);
    // #endif
    bbkj;
 
    lli n, x;cin>>n>>x;
    lli a[n];
    loop(i,0,n-1)
        cin>>a[i];
    vector<lli> prev(x+1, 0), curr(x+1, 0);
    prev[0] = 0;
    loop(i,1,x)
        prev[i] = i%a[0] == 0 ? 1:0;
 
 
    if(n == 1)
        curr = prev;
    // lli ans = 0;
    loop(i,1,n-1){
        curr = vector<lli>(x+1, 0);
        curr[0] = 0;
        loop(j,1,x){
            if(j == a[i])
                curr[j] = 1;
            if(j >= a[i]){
                curr[j] += (prev[j] + curr[j-a[i]])%p;
            }
            else
                curr[j] += prev[j]%p;
        }
        // for(auto i:curr)
        //     cout<<i<<" ";
        // cout<<endl;
 
        prev = curr;
    }    
 
    cout<<curr[x]%p<<endl;
 
    // loop(j,1,x){
    //     loop(i,0,n-1)
    //         cout<<dp[i][j]<<" ";
    //     cout<<endl;
    // }
 
}
