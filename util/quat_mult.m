% Multiply two quaternions    
function result_quat = quat_mult(q1, q2)
    result_quat = zeros(4,1);
    mult_mat = [q1(1) -q1(2) -q1(3) -q1(4);
                q1(2)  q1(1) -q1(4)  q1(3);
                q1(3)  q1(4)  q1(1) -q1(2);
                q1(4) -q1(3)  q1(2)  q1(1);
                ];
            
    result_quat = mult_mat*q2;    
    
end

