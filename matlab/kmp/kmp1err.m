function kmp1err(a)
%KMP1ERR Evaluates error for KMP1.
    b = randb(length(a));
    [x,err] = kmp1(a,b);
    err
    
    hold on;
    plot(sort(lap(a)*x - b),'r');
    refline(0,0);
    refline(0,min(b));
    refline(0,max(b));
    hold off;
end

