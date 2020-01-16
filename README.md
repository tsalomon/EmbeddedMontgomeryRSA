# EmbeddedMontgomeryRSA
An optimized C implementation of the RSA public key encryption using the Montgomery Multiplication algorithm.

## RSA and Modular Exponentiation
The RSA operations for encryption and decryption involve modular exponentiation: X^Y mod M.

  * _C_: Ciphertext
  * _P_: PlainText
  * _e_: Public Exponent
  * _d_: Private Exponent
  * _M_: modulo

  - _C = P<sup>e</sup> mod N_   (encrypt with public key)
  - _P = C<sup>d</sup> mod N_   (decrypt with private key)
  
  The public key is (_e_,_M_) and the private key is (_d_,_M_).

### Generation of Modulus (M)
Choose two prime numbers with half the bit-length of the desired public/private keys.

* _p_: prime
* _q_: prime
* _M = p x q_
  
See [https://www.mobilefish.com/services/rsa_key_generation/rsa_key_generation.php] for more information.
  
### Generation of Public Exponent (e)
The public exponent (_e_) is a small integer. 
Valid choices are 3, 5, 7, 17, 257, or 65537.  
With RSA keys the public exponent is usually 65537.

The requirements for _e_ are: 
* 1 < _e_ < phi(_M_) = (_p_ - 1) * (_q_ - 1)
* _e_ and Euler's phi(_M_) function are coprime.

### Generaion of Private Exponent (d)
The private exponent is generated from primes _p_ and _q_ and the public exponent.
* _d_ = _e_-1 mod (_p_ - 1) * (_q_ - 1)

