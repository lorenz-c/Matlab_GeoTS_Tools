function struct_out = inregion(struct_in, mask, area_id)
    




[XM, YM] = meshgrid(mask.Data.lon, mask.Data.lat);

XM(mask ~= area_id) = NaN;
YM(mask ~= area_id) = NaN;

XM = XM(~isnan(XM));
YM = YM(~isnan(YM));

nr_stns = length(struct_in.Data.stations);


for i = 1:nr_stns
    

