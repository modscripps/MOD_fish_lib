function Profile=epsiProcess_add_shear_heading_to_Profile(Profile)
    Profile.vnav=epsiProcess_cpt_shear_heading(Profile.vnav);
    Profile.shear_heading=interp1(Profile.vnav.dnum, ...
                                  Profile.vnav.shear_heading,...
                                  Profile.dnum);
end