function T = incidencematrix(branch_data)
  % Get the number of buses and transmission lines
  from_bus = branch_data(:,1);
  to_bus = branch_data(:,2);
  num_buses = max([from_bus; to_bus]);
  num_lines = length(from_bus);

  % Initialize the incidence matrix to all zeros
  T = zeros(num_buses, num_lines);

  % Populate the incidence matrix
  for i = 1:num_lines
      T(from_bus(i), i) = 1;
      T(to_bus(i), i) = -1;
  end
end
