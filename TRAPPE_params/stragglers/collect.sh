#!/bin/bash

line=64 

response="$(curl 'http://trappe.oit.umn.edu/scripts/download_parameters.php' \
-H 'Connection: keep-alive' \
-H 'Cache-Control: max-age=0' \
-H 'Origin: http://trappe.oit.umn.edu' \
-H 'Upgrade-Insecure-Requests: 1' \
-H 'DNT: 1' \
-H 'Content-Type: application/x-www-form-urlencoded' \
-H 'User-Agent: Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/86.0.4240.75 Safari/537.36' \
-H 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.9' \
-H 'Referer: http://trappe.oit.umn.edu/' \
-H 'Accept-Language: en-US,en;q=0.9' \
-H 'Cookie: _ga=GA1.2.1049003370.1529510584; _gid=GA1.2.1663596491.1602819527' \
--data-raw "molecule_id=$line" \
--compressed \
--insecure)"
echo "$response" >> trappe_parameters_$line.txt
