.. code-block:: console

    $ pypeit_chk_for_calibs -h
    usage: pypeit_chk_for_calibs [-h] [-s SPECTROGRAPH] [-e EXTENSION]
                                 [--save_setups]
                                 root
    
    Script to check for calibrations [v1]
    
    positional arguments:
      root                  File path+root, e.g. /data/Kast/b
    
    optional arguments:
      -h, --help            show this help message and exit
      -s SPECTROGRAPH, --spectrograph SPECTROGRAPH
                            A valid spectrograph identifier: keck_deimos,
                            keck_lris_blue, keck_lris_red, keck_lris_red_orig,
                            keck_lris_blue_orig, keck_nires, keck_nirspec_low,
                            keck_mosfire, keck_hires_red, keck_kcwi,
                            shane_kast_blue, shane_kast_red, shane_kast_red_ret,
                            tng_dolores, wht_isis_blue, wht_isis_red,
                            vlt_xshooter_uvb, vlt_xshooter_vis, vlt_xshooter_nir,
                            vlt_fors2, gemini_gnirs, gemini_flamingos1,
                            gemini_flamingos2, gemini_gmos_south_ham,
                            gemini_gmos_north_e2v, gemini_gmos_north_ham,
                            magellan_fire, magellan_fire_long, magellan_mage,
                            lbt_mods1r, lbt_mods1b, lbt_mods2r, lbt_mods2b,
                            lbt_luci1, lbt_luci2, mmt_binospec, mmt_mmirs,
                            mdm_osmos_mdm4k, not_alfosc, p200_dbsp_blue,
                            p200_dbsp_red, p200_tspec (default: None)
      -e EXTENSION, --extension EXTENSION
                            File extension; compression indicators (e.g. .gz) not
                            required. (default: .fits)
      --save_setups         If not toggled, remove setup_files/ folder and its
                            files. (default: False)
    