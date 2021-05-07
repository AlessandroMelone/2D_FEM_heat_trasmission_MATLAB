function [] = draw_stabilityDiagram(center,radius, beta)
    
    angle = 0:0.001:2*pi;
    x = center(1) + radius*cos(angle);
    y = center(2) + radius*sin(angle);
    
    Re_beta = real(beta);
    Im_beta = imag(beta);
    
    figure, plot(x,y);
    grid on, xlabel('Re \{B\}'), ylabel('Im \{B\}'), hold on;
    scatter(Re_beta, Im_beta,'x'); title('Stability diagram - Euler expl');

end

