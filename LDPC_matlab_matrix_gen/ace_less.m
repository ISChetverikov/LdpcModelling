function result = ace_less(ace1, ace2)
    big_value = 1e10;
    result = 0;
    ace1(ace1==-1)=big_value;
    ace2(ace2==-1)=big_value;
    temp = ace1-ace2;
    if temp(find(temp,1)) < 0
        result = 1;
    end
end

