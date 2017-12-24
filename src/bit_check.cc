int64_t odd_mask = 0;
int64_t even_mask = 0;

while (int i = 0; i < 32; i++) {
	odd_mask << 2;
	odd_mask |= 2;
	even_mask << 2;
	even_mask |= 1;
}
cout << odd_mask << ", " << even_mask;

inline bool check_hamming_distance(int64_t k, int64_t l){

    // Suppose these two integers converted from two 32mers are k and l:
    x = k ^ l; // to find different bits between k and l 
    
    if(x == 0) return true; // k and l direct match
    
    if(x & (x-1) == 0) return true; // difference between k and l are only one bits
    
    // next, we need to check difference between k and l are only two consecutive bits, higher bit should be in odd position, and lower bit should be in even position
    y = x & odd_mask; // get odd position
    if(y & (y-1) != 0) return false; // only one bit in odd position
    
    z = x & even_mask; // get even position, 
    if(z & (z-1) != 0) return false; // only one bit in even position
    
    if(y == (z << 1)) return true; // y and z should be consecutive
    else return false;

}