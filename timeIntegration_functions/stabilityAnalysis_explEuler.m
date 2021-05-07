function [] = stabilityAnalysis_explEuler(A,Dt)

    beta = eig(A);
    center = [-1/Dt 0]; radius = 1/Dt;
    draw_stabilityDiagram(center, radius, beta);
end

