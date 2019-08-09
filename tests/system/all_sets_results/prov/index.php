
        <?php
        # Taken from:
        # https://stackoverflow.com/questions/3785055/how-can-i-create-a-simple-index-html-file-which-lists-all-files-directories
        $path = ".";
        $dh = opendir($path);
        $i=1;
        while (($file = readdir($dh)) !== false) {
            if($file != "." && $file != ".." && $file != "index.php" && $file != ".htaccess" && $file != "error_log" && $file != "cgi-bin") {
                echo "<a href='$path/$file'>$file</a><br /><br />";
                $i++;
            }
        }
        closedir($dh);
        ?>
        