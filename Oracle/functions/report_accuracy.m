function [] = report_accuracy(name,A,B,tests)

disp(name+":");
fprintf("number of CI tests: %d \n",tests);
[extra_edges,missing_edges,precision,recall,skeleton_F1_score]=...
    learning_errors(A, B);
fprintf('extra edges: %d\nmissing edges: %d\nprecision: %0.2f\nrecall: %0.2f\nF1 score: %0.2f\n\n\n',extra_edges, missing_edges, precision, recall, skeleton_F1_score);
end

