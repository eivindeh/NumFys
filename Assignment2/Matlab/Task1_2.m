clear all
q = sym('q','real');
B = sym('B','real');
m = sym('m','real');
t = sym('t','real');

v0=sym('v0_',[3,1]);
A=[zeros(3) eye(3);
   zeros(3) [0     q*B/m 0;
             -q*B/m 0     0;
             0     0     0]];
v0(2)=0;        
x = expm(A*t)*[0; 0; 0; v0];
