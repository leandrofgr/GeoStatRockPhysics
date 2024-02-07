function [WELLS] = add_welllog_from_cube(WELLS,property,string_name)

for well = 1:length(WELLS)
    WELLS(well).(string_name) = squeeze(property( WELLS(well).i, WELLS(well).j,:));
end

end

