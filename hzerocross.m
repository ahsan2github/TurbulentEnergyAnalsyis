function [n] = hzerocross(sig)
  n = 0;
  for i=1:numel(sig)-1
    if(sig(i)*sig(i+1) < 0)
      n = n +1;
    end
  end
end
