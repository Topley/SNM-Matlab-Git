
for i =1:length(MUFiring)
    decay = [];
    rate = [];
    for j=1:length(MUFiring{i})
    decay(j) = length(MUFiring{i}(j:end));
    rate = diff(MUFiring{i})*fsamp;
    end 
    hold all
    plot(-log(decay(2:end)./rate))
    drawnow
    pause(1)
end 