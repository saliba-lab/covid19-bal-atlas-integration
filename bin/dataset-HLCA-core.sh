#!/usr/bin/env bash

file_out=data/HLCA/core.h5ad
url="https://corpora-data-prod.s3.amazonaws.com/1beed89c-4f17-432f-b58e-d157fdd7ba84/local.h5ad?AWSAccessKeyId=ASIATLYQ5N5X2WMPR2RG&Signature=DE7IK2Zg4f9CHkmP3vKQX0uOlTk%3D&x-amz-security-token=IQoJb3JpZ2luX2VjECoaCXVzLXdlc3QtMiJHMEUCIQCWD2UyYFpa83aY%2FEBndXMyBObxYBZCXgEZZmLC2ddsnQIgSMj7Q65KrPhiBB76BLYriWJOuZ3Zg0RWF9sF1R4HbGkq9AMIo%2F%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FARABGgwyMzE0MjY4NDY1NzUiDB2bFKRdV2su01b3CyrIA6Za4AkI9WGX0A4o5edLs9Gj6CihYOe1yGJwMsj%2BDGuf9wp7jpiGkD4noZHaBShUeoaQysUHzGt9eFv39HKYS0KnpyyaJHl4Je5MYk%2Fb31VrM7t4gUrcNtHyNXSDdW7C7kMQCXD78gcNvpl4e5BZ5rrjIv6ocJ2skf55PG1tfI3aGDwR3Hg17fMcvT4aAUYipWVrBWBLc6%2BHY2QSM1C9VLxj43kuk8QpFezSzeOCpFY6DuIGEy1cr6%2BPGiQRrteSEgYs10PDsP9NKhaReyfd6g7WJKQakj2BwT9mkYAuU2WP4wYs25e5tpF8y2W4tOHDnJEQDBnnp%2BM5yXSLzYIIDFdL5zNZNUDEgoXxvecvs%2F8xHmUjCJX5%2F2jeScFTZXn6Pr5gsVpdHvwmpEzvRC2Pa0r0trESzVjGtgRvToCHu%2BbovPodNRVPOujAhQvJH5%2BnsBe4B3n9%2B1bCSL27gqCH3K1gSYfgVksCvOIi%2BUtDyMen31umaIRatjGXmyHLmUJWCeByfoiggaZHZnkGbEyIylyIj%2FM0OcgXP3LG04Tkk24bkZcwhzzD8ElXxLc10WVB1uThR9jsbEt3P%2BLQfjLvwSkuG4WK1Pk1OTCLz%2BOeBjqlARg%2FEZnLWdDiMiLfJcBs5OX3wN3ZtZPfc%2Fa5m0OBAiXl19fyyDwxpopS367F07sAsfjdQ2f%2FhHpPo8RbZR6ygQODuJTLUwouGmVchlC17Nz3t85%2Byi4oKR%2FqDcssZlqsO7dnmosaWXP2%2FBNAfqqqRTpGySr%2FibzOw1iMnnk9j3EdLAwkfdkMvepoN1pve4PZ%2B33vSmaHAbGPohGSTxAVvCWwUZlkOQ%3D%3D&Expires=1675766252"

dir=$(dirname $file_out)
if [[ -d $dir ]]; then
  echo "Downloading file to $dir"
else
  echo "Creating $dir..."
  mkdir $dir
fi

if [[ -e $file_out ]]; then
  echo "File already exists. Exiting"
  exit 4
fi

echo Download file from $url
curl -o $file_out $url

exit 0
