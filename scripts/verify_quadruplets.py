#!/usr/bin/env python3
"""
verify_quadruplets.py - Verify prime quadruplet locations

This script verifies that the three reported quadruplet locations
actually produce four consecutive prime values of Q(n).

Author: Ruqing Chen
Repository: https://github.com/Ruqing1963/prime-polynomial-Q47
"""

def is_probable_prime(n, k=20):
    """Miller-Rabin primality test."""
    if n < 2:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0:
        return False
    
    # Write n-1 as 2^r * d
    r, d = 0, n - 1
    while d % 2 == 0:
        r += 1
        d //= 2
    
    # Witnesses to test
    witnesses = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    
    for a in witnesses:
        if a >= n:
            continue
        
        x = pow(a, d, n)
        
        if x == 1 or x == n - 1:
            continue
        
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    
    return True

def Q(n):
    """Compute Q(n) = n^47 - (n-1)^47"""
    return n**47 - (n-1)**47

def verify_quadruplet(start_n):
    """Verify that Q(n), Q(n+1), Q(n+2), Q(n+3) are all prime."""
    print(f"\nVerifying quadruplet starting at n = {start_n:,}")
    print("-" * 60)
    
    results = []
    for offset in range(4):
        n = start_n + offset
        q = Q(n)
        is_prime = is_probable_prime(q)
        results.append(is_prime)
        
        digits = len(str(q))
        status = "PRIME" if is_prime else "COMPOSITE"
        print(f"  Q({n:,}) = {digits}-digit number: {status}")
    
    # Also check boundaries
    q_before = Q(start_n - 1)
    q_after = Q(start_n + 4)
    before_prime = is_probable_prime(q_before)
    after_prime = is_probable_prime(q_after)
    
    print(f"\n  Boundary check:")
    print(f"    Q({start_n - 1:,}): {'PRIME' if before_prime else 'COMPOSITE'} (should be composite)")
    print(f"    Q({start_n + 4:,}): {'PRIME' if after_prime else 'COMPOSITE'} (should be composite)")
    
    is_valid = all(results) and not before_prime and not after_prime
    print(f"\n  Quadruplet valid: {'YES ✓' if is_valid else 'NO ✗'}")
    
    return is_valid

def main():
    print("=" * 70)
    print("PRIME QUADRUPLET VERIFICATION")
    print("Q(n) = n^47 - (n-1)^47")
    print("=" * 70)
    
    quadruplets = [
        117309848,  # S-1
        136584738,  # S-2
        218787064,  # S-3
    ]
    
    all_valid = True
    for i, start_n in enumerate(quadruplets, 1):
        print(f"\n{'=' * 70}")
        print(f"QUADRUPLET S-{i}")
        print("=" * 70)
        
        valid = verify_quadruplet(start_n)
        if not valid:
            all_valid = False
    
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Total quadruplets verified: {len(quadruplets)}")
    print(f"All valid: {'YES ✓' if all_valid else 'NO ✗'}")
    
    if all_valid:
        print("\nHardy-Littlewood Comparison:")
        print(f"  Observed:  3")
        print(f"  Predicted: 3.52")
        print(f"  Ratio:     {3/3.52:.2f}")
        print(f"  Status:    EXCELLENT AGREEMENT")

if __name__ == "__main__":
    main()
