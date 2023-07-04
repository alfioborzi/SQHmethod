
function val = delta_1(t1,x,y,h,u_1,u_2,OCPA)
w_1 = h*B_1(t1,x,y,h,u_1,u_2)/C_1(t1,x,y,h,OCPA);

if(sign(w_1) >= 0)
    if(abs(w_1)<10^-12)
        val = 0.5;
    elseif(abs(C_1(t1,x,y,h,OCPA))<10^-12)
        val = 0;
    else
        val = (1/w_1)-(1/(exp(w_1)-1));
    end
else
    if(abs(w_1)<10^-12)
        val = 0.5;
    elseif(abs(C_1(t1,x,y,h,OCPA))<10^-12)
        val = 1;
    else
        val = (1/w_1)-(1/(exp(w_1)-1));
    end
end
end